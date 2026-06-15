"""
Sentaurus Process 入力ファイル → Cprocess 変換器
=================================================

Sentaurus Process (sprocess) の ``.cmd`` ファイルを読み込み、Cprocess の
デッキ文字列（``cp.run_deck`` が解釈する形式）へ翻訳する。対応関係の網羅表は
``docs/sprocess_command_map.md`` を参照。

設計方針:
  * ``pdbSet`` などのパラメータDBは再現しない（警告コメント化）。
  * 移動境界を伴う deposit/etch/oxidation は構造を変えないため翻訳せず警告。
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

# 元素名 → Cprocess の記号
_SPECIES = {
    "boron": "B", "phosphorus": "P", "arsenic": "As", "antimony": "Sb",
    "b": "B", "p": "P", "as": "As", "sb": "Sb",
}

# Sprocess の <unit> → Cprocess 接尾辞
_UNIT = {
    "um": "um", "nm": "nm", "cm": "cm", "m": "m",
    "kev": "keV", "ev": "eV", "mev": "MeV",
    "min": "min", "s": "s", "sec": "s", "hr": "h", "h": "h",
    "c": "C", "k": "K",
    "cm-2": "", "cm-3": "", "cm2": "", "cm3": "",
}


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


def _species(kv: Dict[str, str], bare: List[str], warn: List[str]) -> str:
    """species= / field= / bare 語から Cprocess 記号を得る。"""
    name = kv.get("species") or kv.get("field") or kv.get("name")
    if name is None:
        # implant Boron ... のように bare 語で来る場合
        for b in bare[1:]:
            if b.lower() in _SPECIES:
                name = b
                break
    if name is None:
        return ""
    sym = _SPECIES.get(name.lower())
    if sym is None:
        warn.append(f"unknown species '{name}' (B/P/As/Sb only)")
        return ""
    return sym


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
        sym = _species(kv, bare, warn) or "P"
        conc, _ = _split_unit(kv.get("concentration") or kv.get("conc", "1e15"))
        lines.append(f"init species={sym} conc={conc}")
    return lines


def _cmd_implant(st: _State, bare, kv, warn, out) -> List[str]:
    _emit_mesh(st, out)
    sym = _species(kv, bare, warn)
    if not sym:
        return [f"# [skip] implant: {warn[-1] if warn else 'no species'}"]
    parts = [f"implant species={sym}"]
    if "dose" in kv:
        parts.append(f"dose={_q(kv['dose'])}")
    if "energy" in kv:
        parts.append(f"energy={_q(kv['energy'], 'keV')}")
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
    return [" ".join(parts)]


def _cmd_diffuse(st: _State, bare, kv, warn, out) -> List[str]:
    _emit_mesh(st, out)
    lines = []
    ambient = ("o2" in " ".join(bare).lower() or "h2o" in " ".join(bare).lower()
               or "gas_flow" in kv)
    if ambient:
        lines.append("# [warn] oxidation ambient ignored (Deal-Grove not implemented); "
                     "running inert anneal only")
    parts = ["diffuse"]
    if "time" in kv:
        parts.append(f"time={_q(kv['time'], 'min')}")
    if "temperature" in kv:
        parts.append(f"temp={_q(kv['temperature'], 'C')}")
    elif "temp" in kv:
        parts.append(f"temp={_q(kv['temp'], 'C')}")
    lines.append(" ".join(parts))
    return lines


def _cmd_mask(st: _State, bare, kv, warn) -> List[str]:
    if "left" in kv and "right" in kv:
        st.mask = (st._loc_um(kv["left"]), st._loc_um(kv["right"]))
    return [f"# [info] mask window {st.mask} um (applied to following implant)"]


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


_IGNORE = {"refinebox", "grid", "photo", "strip", "layers", "select",
           "writeplx", "setplxlist", "contact", "pdbset", "pdbget",
           "pdbsetdouble", "pdbdelaydouble", "plot.1d", "plot.2d", "tclsel",
           "graphics", "math"}
_HARD = {"deposit", "etch", "transform", "silicide", "stress", "temp_ramp",
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
        elif name == "mask":
            out += _cmd_mask(st, bare, kv, warn)
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
