"""The dtOO OpenFOAM setup (case/xml/tistos_ru_of.xml) around an existing polyMesh.

    foam_setup.py <case_dir> --dtoo <run>/case [--point n] [--zone-all]

This dtOO build has no OpenFOAM case builder (see case/gen_mesh.py), so the case is
written here from the dtOO XML itself, the same for the dtOO mesh and for ours:
  - system/controlDict, fvSchemes, fvSolution, decomposeParDict and
    constant/transportProperties, turbulenceProperties, MRFProperties: the dictRule
    blocks of the XML, verbatim, with #cV_*# / {ii} substituted;
  - 0/U, p, k, omega, nut: fieldRule + setupRule (per patch) of the XML;
  - constants (#cV_*#) from the templateState of case/machineSave.xml, formulas of
    case/xml/case_constants.xml evaluated;
  - the inlet data (boundaryData_RU_INLET/<point>) is linked next to the case, where
    mappedFieldFixedValue looks for it (../boundaryData_RU_INLET/<point>/boundaryData);
  - --zone-all puts every cell into the cellZone RU of the MRF zone (topoSet).
The patch names must be the dtOO ones (RU_HUB ... RU_PERIOSS).
"""
import argparse
import html
import math
import re
import subprocess
from pathlib import Path

HEADER = ("FoamFile\n{{\n    version     2.0;\n    format      ascii;\n    class       {cls};\n"
          "    object      {obj};\n}}\n\n")
FOAM_BASHRC = "/usr/lib/openfoam/openfoam2606/etc/bashrc"


def constants(case):
    """#cV_*# values: the templateState numbers, formulas from case_constants.xml."""
    st = (case / "machineSave.xml").read_text()
    st = st[st.index('<state label="templateState">'):]
    st = st[:st.index("</state>")]
    val = {}
    for m in re.finditer(r'<constValue[^>]*?value="([^"]*)"[^>]*?label="([^"]*)"', st):
        val[m.group(2)] = m.group(1)

    def ev(name, depth=0):
        v = val[name]
        try:
            return float(v)
        except ValueError:
            pass
        expr = re.sub(r"#(cV_[A-Za-z0-9_]+)#", lambda mm: repr(ev(mm.group(1), depth + 1)), v)
        expr = expr.replace("_pi", repr(math.pi)).replace("^", "**")
        return float(eval(expr, {"__builtins__": {}}))
    return ev


def blocks(text):
    """Top-level `name { ... }` blocks of a dictRule string -> {name: body}."""
    out, i = {}, 0
    while True:
        m = re.compile(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*\{").match(text, i)
        if not m:
            break
        depth, j = 1, m.end()
        while depth:
            depth += {"{": 1, "}": -1}.get(text[j], 0)
            j += 1
        out[m.group(1)] = text[m.end():j - 1]
        i = j
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("case")
    ap.add_argument("--dtoo", required=True, help="the run's dtOO case directory (machineSave.xml, xml/)")
    ap.add_argument("--point", default="n", help="operating point {ii}: n, tl, vl")
    ap.add_argument("--zone-all", action="store_true")
    ap.add_argument("--subdomains", type=int, default=None)
    a = ap.parse_args()
    case, dt = Path(a.case).resolve(), Path(a.dtoo).resolve()
    ev = constants(dt)
    x = (dt / "xml" / "tistos_ru_of.xml").read_text()

    def attr(name):
        i = x.index(name + '="') + len(name) + 2
        return html.unescape(x[i:x.index('"', i)])        # unescape after cutting the attribute

    def subst(s):
        s = s.replace("{ii}", a.point)
        s = re.sub(r"`([^`]*)`", lambda m: repr(float(eval(
            re.sub(r"#(cV_[A-Za-z0-9_]+)#", lambda mm: repr(ev(mm.group(1))), m.group(1))
            .replace("^", "**"), {"__builtins__": {}}))), s)
        return re.sub(r"(-?)#(cV_[A-Za-z0-9_]+)#",
                      lambda m: repr((-1 if m.group(1) else 1) * ev(m.group(2))), s)

    # ---- dictionaries
    for name, body in blocks(attr("dictRule")).items():
        fn = re.search(r'dtOOFileName\s+"([^"]+)"\s*;', body).group(1)
        body = re.sub(r'dtOOFileName\s+"[^"]+"\s*;', "", body)
        body = subst(body)
        if name == "decomposeParDict":                   # no METIS in this OpenFOAM build
            body = re.sub(r"method\s+metis\s*;", "method          scotch;", body)
        if a.subdomains and name == "decomposeParDict":
            body = re.sub(r"numberOfSubdomains\s+\d+", f"numberOfSubdomains {a.subdomains}", body)
        p = case / fn
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(HEADER.format(cls="dictionary", obj=p.name) + body + "\n")

    # ---- fields: internal value + per-patch BCs
    fields = {}
    for m in re.finditer(r":(vol\w+Field)\s*:\s*:?(\w+):\s*:?\[([^\]]*)\]\s*:\s*:([^:]*):", attr("fieldRule")):
        fields[m.group(2)] = (m.group(1), m.group(3), subst(m.group(4)).strip())
    setup = attr("setupRule")
    bc = {f: {} for f in fields}
    rules, depth, start = [], 0, None
    for i, ch in enumerate(setup):                      # top-level { ... } with nesting
        if ch == "{":
            depth += 1
            start = i + 1 if depth == 1 else start
        elif ch == "}":
            depth -= 1
            if depth == 0:
                rules.append(setup[start:i])
    for r in rules:
        head = re.match(r"\s*:(\w+)::([A-Z_,]+):", r)
        if not head:
            continue
        kind, patches = head.group(1), head.group(2).split(",")
        if kind == "ofOpenFOAMCyclicAmiRule":
            for p in patches:
                for f in fields:
                    bc[f][p] = "type cyclicAMI;"
            continue
        for fm in re.finditer(r":(\w+)::(\w+)\((.*?)\):", r[head.end():], re.S):
            f, typ, args = fm.group(1), fm.group(2), subst(fm.group(3)).strip()
            for p in patches:
                bc.setdefault(f, {})[p] = f"type {typ};\n        " + args.replace("\n", "\n        ")
    # patches without a rule for a field (nut on inlet / outlet): calculated, as dtOO
    allp = re.findall(r"^\s{4}(\w+)\s*\n\s{4}\{", (case / "constant" / "polyMesh" / "boundary").read_text(), re.M)
    for f, (cls, dims, iv) in fields.items():
        uni = f"uniform {iv}"
        for p_ in allp:
            bc[f].setdefault(p_, f"type calculated;\n        value {uni};")
        bf = "\n".join(f"    {p}\n    {{\n        {v}\n    }}" for p, v in sorted(bc.get(f, {}).items()))
        (case / "0").mkdir(exist_ok=True)
        (case / "0" / f).write_text(HEADER.format(cls=cls, obj=f) + f"dimensions [{dims}];\n\ninternalField {uni};\n\n"
                                    f"boundaryField\n{{\n{bf}\n}}\n")
    # ---- inlet data next to the case
    link = case.parent / "boundaryData_RU_INLET"
    if not link.exists():
        link.symlink_to(dt / "boundaryData_RU_INLET")
    if a.zone_all:
        (case / "system" / "topoSetDict").write_text(HEADER.format(cls="dictionary", obj="topoSetDict") +
            "actions\n(\n    { name RU; type cellZoneSet; action new; source boxToCell;\n"
            "      box (-1e9 -1e9 -1e9) (1e9 1e9 1e9); }\n);\n")
        subprocess.run(["bash", "-c", f"source {FOAM_BASHRC} >/dev/null 2>&1; cd '{case}'; topoSet > log.topoSet 2>&1"],
                       check=True)
    print(f"[setup] {case}: point {a.point}, omega {ev('cV_omega'):.4f} rad/s, U {ev('cV_U_' + a.point):.3f} m/s, "
          f"fields {sorted(fields)}")


if __name__ == "__main__":
    main()
