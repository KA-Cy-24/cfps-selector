#!/usr/bin/env python3
"""Package the public JS/CSS into an offline-only HTML (never embeds data)."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
css=(ROOT/'repo_additions/static/cfps/model-free/model_free.css').read_text()
js=(ROOT/'repo_additions/static/cfps/model-free/model_free_injector.js').read_text()
if '</script' in js.lower():raise RuntimeError('Unsafe inline script delimiter')
html='''<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>HTGAA CFPS Model-Free Selection — private local viewer</title>
<style>'''+css+'''</style></head><body class="mf-page">
<main class="mf-inner"><p class="mf-muted">Read-only local preview. Choose your private model_free_analysis.json. The file stays in this browser and is not uploaded.</p>
<div id="cfps-model-free-root" class="mf-root"></div></main>
<script>'''+js+'''</script><script>
window.CFPSModelFreeView.mount(document.getElementById('cfps-model-free-root'));
</script></body></html>'''
out=ROOT/'viewer/CFPS_Model_Free_Standalone_Viewer.html';out.write_text(html,encoding='utf-8')
print(out,len(html),'bytes, no private JSON embedded')
if __name__=='__main__':pass
