from pathlib import Path
import re

path = Path(r"C:\Document\Serieux\Travail\Package\git\susieR\vignettes\susie_slide_prior.tex")
text = path.read_text(encoding="utf-8")
labels = re.findall(r"\\label\{([^}]+)\}", text)
refs = re.findall(r"\\(?:eqref|ref)\{([^}]+)\}", text)
assert len(labels) == len(set(labels)), "Duplicate equation/section labels"
assert not set(refs) - set(labels), ("Missing references", set(refs)-set(labels))
keys = set(re.findall(r"\\bibitem\{([^}]+)\}", text))
cited = set(k for group in re.findall(r"\\cite\{([^}]+)\}", text) for k in group.split(","))
assert cited <= keys, ("Missing citations", cited-keys)
clean = re.sub(r"\\begin\{verbatim\}.*?\\end\{verbatim\}", "", text, flags=re.S)
clean = re.sub(r"(?<!\\)%[^\n]*", "", clean)
stack = []
for match in re.finditer(r"\\(begin|end)\{([^}]+)\}", clean):
    action, env = match.groups()
    if action == "begin":
        stack.append(env)
    else:
        assert stack and stack.pop() == env, f"Unmatched environment {env}"
assert not stack, stack
brace_level = 0
for match in re.finditer(r"(?<!\\)[{}]", clean):
    brace_level += 1 if match.group() == "{" else -1
    assert brace_level >= 0, "Closing brace without opening brace"
assert brace_level == 0, f"Brace balance {brace_level}"
assert len(re.findall(r"(?<!\\)\$", clean)) % 2 == 0, "Unmatched inline math delimiter"
assert "^{,2}" not in text
assert "initialization; if they are fixed permanently" not in text
assert text.count(r"\end{document}") == 1
print(f"Static checks passed: {len(labels)} unique labels, {len(refs)} references, {len(keys)} bibliography entries.")
print("LaTeX environments, braces and inline math delimiters balance.")
print("These checks do not establish successful typesetting or validate page layout.")
