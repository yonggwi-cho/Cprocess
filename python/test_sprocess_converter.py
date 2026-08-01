"""Golden tests for cprocess.sprocess (W-1 converter overhaul).

Run standalone: python3 python/test_sprocess_converter.py
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from cprocess.sprocess import translate  # noqa: E402

_FAILS = 0


def check(name, cond, detail=""):
    global _FAILS
    status = "PASS" if cond else "FAIL"
    print(f"[{status}] {name}" + (f" - {detail}" if detail and not cond else ""))
    if not cond:
        _FAILS += 1


def assert_no_unsupported(deck: str, label: str):
    lines = [l for l in deck.splitlines() if "[unsupported]" in l]
    check(f"{label}: no [unsupported] lines", not lines, "; ".join(lines))


def assert_deck_syntax(deck: str, label: str):
    """Every non-comment, non-empty line's key=value tokens parse as
    key=value (matches src/deck.cpp's simple tokenizer contract)."""
    ok = True
    bad = []
    for line in deck.splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        toks = line.split()
        for t in toks[1:]:
            if "=" in t:
                k, v = t.split("=", 1)
                if not k or not v:
                    ok = False
                    bad.append(line)
    check(f"{label}: deck syntax valid", ok, "; ".join(bad))


# ---------------------------------------------------------------------------
# 1. LOCOS + well + S/D flow (includes BF2 decomposition)
# ---------------------------------------------------------------------------
LOCOS_RECIPE = """
line x location=0<um> spacing=0.01<um>
line x location=1<um> spacing=0.01<um>
line y location=0<um>
line y location=1<um>
region silicon xlo=0 xhi=1
init concentration=1e15 field=Boron

implant Boron dose=5e12<cm-2> energy=40<keV>

mask name=well left=0<um> right=0.4<um>
implant Phosphorus dose=2e13<cm-2> energy=150<keV>

diffuse temperature=1000<C> time=200<min> O2

mask name=sd left=0<um> right=0.35<um>
implant BF2 dose=4e15<cm-2> energy=15<keV>

diffuse temperature=950<C> time=10<min>
struct tdr=nmos.tdr
exit
"""


def test_locos_well_sd():
    deck = translate(LOCOS_RECIPE)
    assert_no_unsupported(deck, "LOCOS+well+S/D")
    assert_deck_syntax(deck, "LOCOS+well+S/D")
    check("LOCOS: oxidize present", "oxidize" in deck, deck)
    check("LOCOS: BF2 decomposed to B", "implant species=B" in deck, deck)
    check("LOCOS: BF2 energy scaled down from 15keV",
          re.search(r"implant species=B dose=4e15 energy=3\.3\d*keV", deck) is not None,
          deck)


# ---------------------------------------------------------------------------
# 2. diffuse + O2 (dry oxidation)
# ---------------------------------------------------------------------------
DIFFUSE_O2_RECIPE = """
line x location=0.5<um>
init concentration=1e15 field=Boron
diffuse temperature=1000<C> time=30<min> O2
"""


def test_diffuse_o2_dry():
    deck = translate(DIFFUSE_O2_RECIPE)
    assert_no_unsupported(deck, "diffuse+O2")
    assert_deck_syntax(deck, "diffuse+O2")
    check("dry ox: oxidize line present",
          any(l.startswith("oxidize ") for l in deck.splitlines()), deck)
    check("dry ox: ambient=dry", "ambient=dry" in deck, deck)
    check("dry ox: no inert-anneal degradation warning",
          "inert anneal" not in deck, deck)
    check("dry ox: time/temp carried over",
          "time=30min" in deck and "temp=1000C" in deck, deck)


# ---------------------------------------------------------------------------
# 3. diffuse + H2O (wet oxidation)
# ---------------------------------------------------------------------------
DIFFUSE_H2O_RECIPE = """
line x location=0.5<um>
init concentration=1e15 field=Boron
diffuse temperature=1050<C> time=45<min> H2O
"""


def test_diffuse_h2o_wet():
    deck = translate(DIFFUSE_H2O_RECIPE)
    assert_no_unsupported(deck, "diffuse+H2O")
    assert_deck_syntax(deck, "diffuse+H2O")
    check("wet ox: oxidize line present",
          any(l.startswith("oxidize ") for l in deck.splitlines()), deck)
    check("wet ox: ambient=wet", "ambient=wet" in deck, deck)


# ---------------------------------------------------------------------------
# 4. deposit + etch pair (rate*time thickness conversion)
# ---------------------------------------------------------------------------
DEPOSIT_ETCH_RECIPE = """
line x location=0.5<um>
init concentration=1e15 field=Boron
deposit material=oxide thickness=0.1<um>
deposit material=nitride rate=0.02<um/min> time=10<min>
etch material=nitride thickness=0.05<um>
etch material=oxide rate=0.01<um/min> time=20<min>
"""


def test_deposit_etch():
    deck = translate(DEPOSIT_ETCH_RECIPE)
    assert_no_unsupported(deck, "deposit+etch")
    assert_deck_syntax(deck, "deposit+etch")
    check("deposit: thickness passthrough",
          "deposit material=oxide thickness=0.1um" in deck, deck)
    check("deposit: rate*time = 0.2um",
          "deposit material=nitride thickness=0.2um" in deck, deck)
    check("etch: thickness passthrough",
          "etch depth=0.05um material=nitride" in deck, deck)
    check("etch: rate*time = 0.2um",
          "etch depth=0.2um material=oxide" in deck, deck)


# ---------------------------------------------------------------------------
# 5. pdbSet override (known key translated, unknown key warned)
# ---------------------------------------------------------------------------
PDBSET_RECIPE = """
line x location=0.5<um>
pdbSet Silicon Boron sper.actfactor 5.0
pdbSet Silicon Boron some.made.up.key 1.0
"""


def test_pdbset():
    deck = translate(PDBSET_RECIPE)
    assert_deck_syntax(deck, "pdbSet")
    check("pdbSet: known key translated",
          "pdbset key=sper.act_factor value=5.0" in deck, deck)
    check("pdbSet: unknown key warned, not silently dropped",
          any("[warn]" in l and "some.made.up.key" in l for l in deck.splitlines()),
          deck)
    check("pdbSet: unknown key does not appear as a pdbset command",
          "key=some.made.up.key" not in deck, deck)


def main():
    test_locos_well_sd()
    test_diffuse_o2_dry()
    test_diffuse_h2o_wet()
    test_deposit_etch()
    test_pdbset()
    print()
    if _FAILS:
        print(f"{_FAILS} FAILURE(S)")
        return 1
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
