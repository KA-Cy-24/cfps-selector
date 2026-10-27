#!/usr/bin/env python3
"""Create the source-only standalone viewer (no participant data embedded)."""
from pathlib import Path
r=Path(__file__).resolve().parents[1]
css=(r/'repo_additions/static/cfps/decision/decision.css').read_text()
js=(r/'repo_additions/static/cfps/decision/decision.js').read_text()
if '</script>' in js.lower():raise ValueError('Unsafe inline script closing tag')
html='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>HTGAA · CFPS Plate Decision Lab</title><style>'+css+'</style></head><body style="margin:0"><div id="cfps-root"></div><script>'+js+'</script><script>window.CFPSDecisionLab.mount(document.getElementById("cfps-root"));</script></body></html>'
p=r/'viewer/CFPS_Decision_Lab_Standalone.html'
p.write_text(html,encoding='utf-8')
print(p,'bytes',p.stat().st_size)
