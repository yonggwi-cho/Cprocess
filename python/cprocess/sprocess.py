"""
Sentaurus Process 入力ファイル → Cprocess 変換器
=================================================

Sentaurus Process (sprocess) の ``.cmd`` ファイルを読み込み、Cprocess の
デッキ文字列（``cp.run_deck`` が解釈する形式）へ翻訳する。対応関係の網羅表は
``docs/sprocess_command_map.md`` を参照。

設計方針 (W-1 で全面更新):
  * ``diffuse ... O2|H2O`` は Deal-Grove 酸化 (``oxidize``/``oxidize2d``) へ
    翻訳する（雰囲気を無視した不活性アニールへの劣化はしない）。
  * ``deposit``/``etch``/``silicide``/``photo``/``strip``/``temp_ramp`` は
    現行デッキ実装 (``src/deck.cpp``) へ直接翻訳する。
  * ``pdbSet`` は既知キー対訳表 (``_PDB_KEY_MAP``) で ``pdbset`` へ翻訳し、
    対応外キーは警告コメントを残す（黙って捨てない）。
  * ドーパント種は B/P/As/Sb/In/C/F/Ge の8種 + BF2（B への質量分解近似）。
  * Sprocess の x（深さ方向）を Cprocess の z（深さ）に対応付ける。
  * 単位は ``<unit>`` 角括弧を剥がして Cprocess 接尾辞へ変換する。

使い方::

    from cprocess.sprocess import translate_file, run_sprocess
    deck = translate_file("nmos.cmd")     # 翻訳のみ
    st   = run_sprocess("nmos.cmd")       # 翻訳して実行

CLI::

    python -m cprocess.sprocess nmos.cmd [-o out.deck] [--run]
"""

from __future__ import annotations

import re
import sys
from typing import Dict, List, Tuple

__all__ = ["translate", "translate_file", "run_sprocess"]

# 元素名/記号 → Cprocess の記号（8種: B/P/As/Sb/In/C/F/Ge）
_SPECIES = {
    "boron": "B", "phosphorus": "P", "arsenic": "As", "antimony": "Sb",
    "indium": "In", "carbon": "C", "fluorine": "F", "germanium": "Ge",
    "b": "B", "p": "P", "as": "As", "sb": "Sb",
    "in": "In", "c": "C", "f": "F", "ge": "Ge",
}

# BF2+ 分子イオン注入の質量ベース分解近似:
#   BF2+ は打ち込み後 B と F に解離するため、Cprocess には無い分子種を
#   実効的に "B原子1個 (ドーズはそのまま)" として扱う。エネルギーは
#   E_B = E_BF2 * (M_B / M_BF2) の質量比スケーリング（SRIM 等で使われる
#   一次近似と同じ考え方）で換算する。電荷分配・分子解離ダイナミクスは
#   無視した近似であり、実プロファイルとの厳密な一致は保証しない。
_M_BF2 = 48.6   # amu (B=10.8 + 2*F=19.0, natural-abundance-ish approx)
_M_B = 11.0     # amu
_BF2_ENERGY_SCALE = _M_B / _M_BF2   # ~= 0.226

# Sprocess の <unit> → Cprocess 接尾辞
_UNIT = {
    "um": "um", "nm": "nm", "cm": "cm", "m": "m",
    "kev": "keV", "ev": "eV", "mev": "MeV",
    "min": "min", "s": "s", "sec": "s", "hr": "h", "h": "h",
    "c": "C", "k": "K",
    "cm-2": "", "cm-3": "", "cm2": "", "cm3": "",
}

# pdbSet の Sentaurus 慣用パラメータ名 (小文字, 記号除去して比較) →
# Cprocess ParamDB キー。値は format 文字列で <mat>/<sym> をプレースホルダ
# として埋め込む場合がある。key 側は英数字のみに正規化して比較する
# (大文字小文字・"_"・"." の有無を無視)。
def _norm_key(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", s.lower())


_PDB_KEY_MAP_STATIC = {
    _norm_key("ted.frenkelsurvival"): "ted.frenkel_survival",
    _norm_key("frenkelsurvival"): "ted.frenkel_survival",
    _norm_key("ted.kci"): "ted.k_ci",
    _norm_key("interstitial.kci"): "ted.k_ci",
    _norm_key("oed.theta"): "oed.theta",
    _norm_key("oxidant.theta"): "oed.theta",
    _norm_key("oed.psicap"): "oed.psi_cap",
    _norm_key("oed.psi_cap"): "oed.psi_cap",
    _norm_key("ox2d.cgas"): "ox2d.cgas",
    _norm_key("oxidant.cgas"): "ox2d.cgas",
    _norm_key("ox2d.nitrideleak"): "ox2d.nitride_leak",
    _norm_key("ox2d.seedox"): "ox2d.seed_ox_um",
    _norm_key("sige.couple"): "sige.couple",
    _norm_key("sige.degcoef"): "sige.dEg_coef",
    _norm_key("sige.bandgapnarrowing"): "sige.dEg_coef",
    _norm_key("sige.eps0coef"): "sige.eps0_coef",
    _norm_key("sper.actfactor"): "sper.act_factor",
    _norm_key("sper.amorphdensity"): "sper.amorph_density",
    _norm_key("sper.ea"): "sper.ea",
    _norm_key("sper.v0"): "sper.v0",
    _norm_key("stress.couple"): "stress.couple",
    _norm_key("stress.vr"): "stress.vr",
}

# <mat>/<sym> プレースホルダを含むキー（正規表現で先頭一致）
_PDB_KEY_MAP_PATTERNS = [
    (re.compile(r"^stressvact(?P<mat>[a-z0-9]+)$"), "stress.vact.{mat}"),
    (re.compile(r"^mech[e](?P<mat>[a-z0-9]+)$"), "mech.E.{mat}"),
    (re.compile(r"^youngsmodulus(?P<mat>[a-z0-9]+)$"), "mech.E.{mat}"),
    (re.compile(r"^mechnu(?P<mat>[a-z0-9]+)$"), "mech.nu.{mat}"),
    (re.compile(r"^poissonratio(?P<mat>[a-z0-9]+)$"), "mech.nu.{mat}"),
    (re.compile(r"^mechalpha(?P<mat>[a-z0-9]+)$"), "mech.alpha.{mat}"),
    (re.compile(r"^mechsigma0(?P<mat>[a-z0-9]+)$"), "mech.sigma0.{mat}"),
    (re.compile(r"^intrinsicstress(?P<mat>[a-z0-9]+)$"), "mech.sigma0.{mat}"),
    (re.compile(r"^mechtau(?P<mat>[a-z0-9]+)$"), "mech.tau.{mat}"),
]


def _pdb_translate(name: str) -> str | None:
    """Sentaurus PDB パラメータ名 → Cprocess ParamDB キー。無ければ None。"""
    nk = _norm_key(name)
    if nk in _PDB_KEY_MAP_STATIC:
        return _PDB_KEY_MAP_STATIC[nk]
    # <Sym>.fi (species-specific)
    m = re.match(r"^([a-z]+)fi$", nk)
    if m:
        sym = _SPECIES.get(m.group(1))
        if sym:
            return f"{sym}.fi"
    for pat, tmpl in _PDB_KEY_MAP_PATTERNS:
        mm = pat.match(nk)
        if mm:
            return tmpl.format(mat=mm.group("mat"))
    return None


def _strip_comments(line: str) -> str:
    """Tcl/Sprocess コメント（先頭または ; の後の #）を除去。"""
    # 行頭 # と末尾 ;# / # コメント
    out, in_str = [], False
    i = 0
    while i < len(line):
        c = line[i]
        if c == '"':
            in_str = not in_str
        if c == "#" and not in_str:
            break
        out.append(c)
        i += 1
    return "".join(out)


def _split_unit(value: str) -> Tuple[str, str]:
    """``20<keV>`` → ('20', 'keV')、``0.5`` → ('0.5', '')。"""
    m = re.match(r"^([^<]*)(?:<([^>]*)>)?$", value.strip())
    if not m:
        return value, ""
    num, unit = m.group(1), (m.group(2) or "")
    suf = _UNIT.get(unit.lower().strip(), "")
    return num, suf


def _q(value: str, default_unit: str = "") -> str:
    """Sprocess の値を Cprocess の数量文字列へ（単位付与）。"""
    num, unit = _split_unit(value)
    if not unit:
        unit = default_unit
    return f"{num}{unit}"


def _f(value: str) -> float:
    num, _ = _split_unit(value)
    return float(num)


def _parse_kv(tokens: List[str]) -> Tuple[List[str], Dict[str, str]]:
    """トークン列を bare 語と key=value に分離。"""
    bare: List[str] = []
    kv: Dict[str, str] = {}
    for t in tokens:
        if "=" in t and not t.startswith("="):
            k, v = t.split("=", 1)
            kv[k.lower()] = v
        else:
            bare.append(t)
    return bare, kv


def _species(kv: Dict[str, str], bare: List[str], warn: List[str]) -> Tuple[str, float | None]:
    """species= / field= / bare 語から Cprocess 記号を得る。

    戻り値: (記号, エネルギースケール or None)。BF2 の場合のみ
    エネルギースケール (~0.226) を返す。それ以外は None（スケールなし）。
    """
    name = kv.get("species") or kv.get("field") or kv.get("name")
    if name is None:
        # implant Boron ... のように bare 語で来る場合
        for b in bare[1:]:
            if b.lower() in _SPECIES or b.lower() == "bf2":
                name = b
                break
    if name is None:
        return "", None
    low = name.lower()
    if low == "bf2":
        warn.append(
            "BF2 decomposed to B (mass-scaled energy approx, "
            f"factor={_BF2_ENERGY_SCALE:.4g}); dose taken as B atom count, "
            "dissociation/charge effects ignored"
        )
        return "B", _BF2_ENERGY_SCALE
    sym = _SPECIES.get(low)
    if sym is None:
        warn.append(f"unknown species '{name}' (B/P/As/Sb/In/C/F/Ge/BF2 only)")
        return "", None
    return sym, None


class _State:
    """変換中に保持する状態（line 範囲・mask 窓など）。"""

    def __init__(self) -> None:
        self.x_locs: List[float] = []   # Sprocess x(深さ) [um]
        self.y_locs: List[float] = []
        self.z_locs: List[float] = []
        self.materials: List[str] = []
        self.mask: Tuple[float, float] | None = None   # (left,right) [um]
        self.tcl_vars: Dict[str, str] = {}
        self.mesh_emitted = False

    def _loc_um(self, value: str) -> float:
        num, unit = _split_unit(value)
        v = float(num)
        if unit == "nm":
            v *= 1e-3
        elif unit == "cm":
            v *= 1e4
        return v   # um


# ---- 各コマンドの翻訳 -------------------------------------------------------

def _emit_mesh(st: _State, out: List[str]) -> None:
    """line で集めた範囲から mesh box を一度だけ出力。"""
    if st.mesh_emitted:
        return
    xmax = max(st.x_locs) if st.x_locs else 1.0      # 深さ[um] → Cprocess z
    ymax = max(st.y_locs) if st.y_locs else 0.1
    zmax_lateral = max(st.z_locs) if st.z_locs else 0.1
    # Sprocess x(深さ) → Cprocess z(深さ)。横方向は薄い箱で近似。
    nz = max(20, int(round(xmax / 0.01)))            # ~10nm 刻み
    out.append(
        f"mesh box xmax={ymax:g}um ymax={zmax_lateral:g}um zmax={xmax:g}um "
        f"nx=8 ny=8 nz={min(nz, 200)}"
    )
    st.mesh_emitted = True


def _cmd_line(st: _State, bare, kv) -> List[str]:
    axis = bare[1].lower() if len(bare) > 1 else "x"
    if "location" in kv:
        loc = st._loc_um(kv["location"])
        {"x": st.x_locs, "y": st.y_locs, "z": st.z_locs}.get(axis, st.x_locs).append(loc)
    return []   # mesh は最初のプロセスstep直前にまとめて出力


def _cmd_region(st: _State, bare, kv, warn) -> List[str]:
    mat = (kv.get("material") or (bare[1] if len(bare) > 1 else "silicon")).lower()
    st.materials.append(mat)
    return []   # region は mesh 出力後に material 指定へ（簡略: 全体 silicon 既定）


def _cmd_init(st: _State, bare, kv, warn, out) -> List[str]:
    _emit_mesh(st, out)
    lines = []
    if "concentration" in kv or "conc" in kv:
        sym, _scale = _species(kv, bare, warn)
        sym = sym or "P"
        conc, _ = _split_unit(kv.get("concentration") or kv.get("conc", "1e15"))
        lines.append(f"init species={sym} conc={conc}")
    return lines


def _cmd_implant(st: _State, bare, kv, warn, out) -> List[str]:
    _emit_mesh(st, out)
    sym, escale = _species(kv, bare, warn)
    lines = [f"# [info] {w}" for w in warn]
    if not sym:
        lines.append(f"# [skip] implant: {warn[-1] if warn else 'no species'}")
        return lines
    parts = [f"implant species={sym}"]
    if "dose" in kv:
        parts.append(f"dose={_q(kv['dose'])}")
    if "energy" in kv:
        e_num, e_unit = _split_unit(kv["energy"])
        e_val = float(e_num)
        if escale is not None:
            e_val *= escale
        if not e_unit:
            e_unit = "keV"
        parts.append(f"energy={e_val:g}{e_unit}")
    if "tilt" in kv:
        parts.append(f"tilt={_split_unit(kv['tilt'])[0]}")
    if "rotation" in kv:
        parts.append(f"rotation={_split_unit(kv['rotation'])[0]}")
    # MC 判定
    bl = " ".join(bare).lower()
    if "sentaurus.mc" in bl or "crystaltrim" in bl or kv.get("type", "").lower() == "mc":
        parts.append("method=mc")
        if "particles" in kv:
            parts.append(f"ions={_split_unit(kv['particles'])[0]}")
    # mask 窓
    if st.mask is not None:
        l, r = st.mask
        parts.append(f"x1={l:g}um x2={r:g}um y1=0um y2=1um")
    lines.append(" ".join(parts))
    return lines


def _detect_oxidation(bare, kv) -> Tuple[bool, bool]:
    """(is_oxidation, is_wet) を返す。"""
    bl = " ".join(bare).lower()
    gas_flow = kv.get("gas_flow", "").lower()
    wet = "h2o" in bl or "steam" in bl or "wet" in bl or "h2o" in gas_flow or "steam" in gas_flow
    dry = "o2" in bl or "dry" in bl or "o2" in gas_flow
    return (wet or dry), wet


def _cmd_diffuse(st: _State, bare, kv, warn, out) -> List[str]:
    _emit_mesh(st, out)
    lines = []
    is_ox, is_wet = _detect_oxidation(bare, kv)
    if is_ox:
        parts = ["oxidize2d" if kv.get("lateral", "").lower() in ("1", "true", "on", "yes")
                 else "oxidize"]
        if "time" in kv:
            parts.append(f"time={_q(kv['time'], 'min')}")
        if "temperature" in kv:
            parts.append(f"temp={_q(kv['temperature'], 'C')}")
        elif "temp" in kv:
            parts.append(f"temp={_q(kv['temp'], 'C')}")
        parts.append(f"ambient={'wet' if is_wet else 'dry'}")
        lines.append(" ".join(parts))
        return lines
    parts = ["diffuse"]
    if "time" in kv:
        parts.append(f"time={_q(kv['time'], 'min')}")
    if "temperature" in kv:
        parts.append(f"temp={_q(kv['temperature'], 'C')}")
    elif "temp" in kv:
        parts.append(f"temp={_q(kv['temp'], 'C')}")
    lines.append(" ".join(parts))
    return lines


def _cmd_temp_ramp(st: _State, bare, kv, warn, out) -> List[str]:
    """temp_ramp name= time= temperature= (複数行の集合) は Sentaurus では
    通常 name= で定義した後 diffuse temp_ramp=<name> で参照するパターンが
    多いが、単純化のため単発 temp_ramp 行に breakpoints= 形式
    (`0:900,10<min>:1050<C>,...`) が来る簡易記法をサポートする。breakpoints
    が無い場合は変換不能として警告する。
    """
    _emit_mesh(st, out)
    bps = kv.get("breakpoints") or kv.get("points")
    if not bps:
        return [f"# [warn] temp_ramp: no breakpoints=/points= given, cannot translate: "
                f"{' '.join(bare)}"]
    # breakpoints="0:900C,10min:1050C,20min:900C"
    segs = []
    for tok in bps.replace(";", ",").split(","):
        tok = tok.strip()
        if not tok:
            continue
        t_s, temp_s = tok.split(":", 1)
        t_num, t_unit = _split_unit(t_s)
        if not t_unit:
            t_unit = "min"
        temp_num, temp_unit = _split_unit(temp_s)
        if not temp_unit:
            temp_unit = "C"
        segs.append(f"{t_num}{t_unit}:{temp_num}{temp_unit}")
    return [f"diffuse ramp={','.join(segs)}"]


def _cmd_mask(st: _State, bare, kv, warn) -> List[str]:
    if "left" in kv and "right" in kv:
        st.mask = (st._loc_um(kv["left"]), st._loc_um(kv["right"]))
    return [f"# [info] mask window {st.mask} um (applied to following implant)"]


def _cmd_photo(st: _State, bare, kv, warn, out) -> List[str]:
    _emit_mesh(st, out)
    thickness = kv.get("thickness", "0.4<um>")
    return [f"photo resist={_q(thickness, 'um')}"]


def _cmd_strip(st: _State, bare, kv, warn, out) -> List[str]:
    _emit_mesh(st, out)
    return ["strip"]


def _rate_time_thickness(kv: Dict[str, str], warn: List[str]) -> str | None:
    """thickness= があればそのまま、rate=+time= なら thickness=rate*time
    (um/min * min = um) を計算して um 単位の量文字列を返す。"""
    if "thickness" in kv:
        return _q(kv["thickness"], "um")
    if "rate" in kv and "time" in kv:
        rate_num, rate_unit = _split_unit(kv["rate"])
        time_num, time_unit = _split_unit(kv["time"])
        # rate は um/min 前提（Sentaurus 既定に合わせる）。time は min 既定。
        rate_val = float(rate_num)
        time_val = float(time_num)
        if time_unit == "h":
            time_val *= 60.0
        elif time_unit == "s":
            time_val /= 60.0
        thickness_um = rate_val * time_val
        return f"{thickness_um:g}um"
    warn.append("deposit/etch: need thickness= or rate=+time=")
    return None


def _cmd_deposit(st: _State, bare, kv, warn, out) -> List[str]:
    _emit_mesh(st, out)
    mat = kv.get("material", "oxide")
    thick = _rate_time_thickness(kv, warn)
    lines = [f"# [info] {w}" for w in warn]
    if thick is None:
        lines.append(f"# [unsupported] deposit: {' '.join(bare)}")
        return lines
    lines.append(f"deposit material={mat} thickness={thick}")
    return lines


def _cmd_etch(st: _State, bare, kv, warn, out) -> List[str]:
    _emit_mesh(st, out)
    mat = kv.get("material", "")
    thick = _rate_time_thickness(kv, warn)
    lines = [f"# [info] {w}" for w in warn]
    if thick is None:
        lines.append(f"# [unsupported] etch: {' '.join(bare)}")
        return lines
    parts = [f"etch depth={thick}"]
    if mat:
        parts.append(f"material={mat}")
    lines.append(" ".join(parts))
    return lines


def _cmd_silicide(st: _State, bare, kv, warn, out) -> List[str]:
    _emit_mesh(st, out)
    metal = kv.get("metal", "nickel")
    parts = [f"silicide metal={metal}"]
    if "time" in kv:
        parts.append(f"time={_q(kv['time'], 'min')}")
    if "temperature" in kv:
        parts.append(f"temp={_q(kv['temperature'], 'C')}")
    elif "temp" in kv:
        parts.append(f"temp={_q(kv['temp'], 'C')}")
    return [" ".join(parts)]


def _cmd_pdbset(st: _State, bare, kv, warn) -> List[str]:
    # pdbSet <Reg> <Sp> <Param> <val>  もしくは key= value= 形式の両対応
    name = kv.get("param") or kv.get("key")
    value = kv.get("value")
    if name is None or value is None:
        # bare 語形式: pdbSet <Reg> <Sp> <Param> <val>
        if len(bare) >= 4:
            name, value = bare[-2], bare[-1]
        elif len(bare) >= 3:
            name, value = bare[-2], bare[-1]
    if name is None or value is None:
        return [f"# [warn] pdbSet: cannot parse '{' '.join(bare)}' — dropped"]
    key = _pdb_translate(name)
    if key is None:
        return [f"# [warn] pdbSet '{name}' has no known Cprocess ParamDB "
                f"equivalent — dropped"]
    val_num, _ = _split_unit(value)
    return [f"pdbset key={key} value={val_num}"]


def _cmd_struct(st: _State, bare, kv, warn) -> List[str]:
    f = kv.get("tdr") or kv.get("gmsh") or kv.get("smesh") or "out.tdr"
    f = re.sub(r"\.(tdr|smesh)$", ".vtu", f)
    if not f.endswith(".vtu"):
        f += ".vtu"
    return [f"save file={f}"]


def _cmd_set(st: _State, bare, kv, warn) -> List[str]:
    # set var value
    if len(bare) >= 3:
        st.tcl_vars[bare[1]] = bare[2]
    return []


def _unsupported(name: str, raw: str) -> List[str]:
    return [f"# [unsupported] {name}: {raw.strip()}"]


_IGNORE = {"refinebox", "grid", "layers", "select",
           "writeplx", "setplxlist", "contact", "pdbget",
           "pdbsetdouble", "pdbdelaydouble", "plot.1d", "plot.2d", "tclsel",
           "graphics", "math"}
_HARD = {"transform",
         "if", "else", "endif", "foreach", "while", "expr"}


def translate(text: str) -> str:
    """Sprocess コマンドテキストを Cprocess デッキ文字列へ翻訳。"""
    st = _State()
    out: List[str] = ["# auto-translated from Sentaurus Process by cprocess.sprocess"]
    for raw in text.splitlines():
        line = _strip_comments(raw).strip()
        if not line:
            continue
        # $var 展開（単純置換）
        for k, v in st.tcl_vars.items():
            line = line.replace(f"${k}", v).replace(f"${{{k}}}", v)
        tokens = line.split()
        if not tokens:
            continue
        name = tokens[0].lower()
        bare, kv = _parse_kv(tokens)
        warn: List[str] = []

        if name == "line":
            out += _cmd_line(st, bare, kv)
        elif name == "region":
            out += _cmd_region(st, bare, kv, warn)
        elif name == "init":
            out += _cmd_init(st, bare, kv, warn, out)
        elif name == "implant":
            out += _cmd_implant(st, bare, kv, warn, out)
        elif name in ("diffuse", "anneal"):
            out += _cmd_diffuse(st, bare, kv, warn, out)
        elif name == "temp_ramp":
            out += _cmd_temp_ramp(st, bare, kv, warn, out)
        elif name == "mask":
            out += _cmd_mask(st, bare, kv, warn)
        elif name == "photo":
            out += _cmd_photo(st, bare, kv, warn, out)
        elif name == "strip":
            out += _cmd_strip(st, bare, kv, warn, out)
        elif name == "deposit":
            out += _cmd_deposit(st, bare, kv, warn, out)
        elif name == "etch":
            out += _cmd_etch(st, bare, kv, warn, out)
        elif name == "silicide":
            out += _cmd_silicide(st, bare, kv, warn, out)
        elif name == "pdbset":
            out += _cmd_pdbset(st, bare, kv, warn)
        elif name == "struct":
            out += _cmd_struct(st, bare, kv, warn)
        elif name == "set":
            out += _cmd_set(st, bare, kv, warn)
        elif name in ("exit", "quit"):
            out.append("stop")
        elif name in _IGNORE:
            out.append(f"# [ignored] {raw.strip()}")
        elif name in _HARD:
            out += _unsupported(name, raw)
        else:
            out += _unsupported(name, raw)
    return "\n".join(out) + "\n"


def translate_file(path: str) -> str:
    """Sprocess ``.cmd`` ファイルを翻訳して Cprocess デッキ文字列を返す。"""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return translate(fh.read())


def run_sprocess(path: str, state=None):
    """Sprocess ファイルを翻訳して即実行し、SimState を返す。"""
    import cprocess as cp
    deck = translate_file(path)
    st = state if state is not None else cp.SimState()
    cp.run_deck(deck, st)
    return st


def _main(argv: List[str]) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Sentaurus Process → Cprocess 変換器")
    ap.add_argument("input", help="Sprocess .cmd ファイル")
    ap.add_argument("-o", "--output", help="デッキ出力先（既定: 標準出力）")
    ap.add_argument("--run", action="store_true", help="変換して実行する")
    args = ap.parse_args(argv)

    deck = translate_file(args.input)
    if args.run:
        import cprocess as cp
        st = cp.SimState()
        print(cp.run_deck(deck, st))
        return 0
    if args.output:
        with open(args.output, "w", encoding="utf-8") as fh:
            fh.write(deck)
    else:
        sys.stdout.write(deck)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))
