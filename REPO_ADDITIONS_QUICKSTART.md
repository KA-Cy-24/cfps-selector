**For the integrated plate decision framework, follow `COMMIT_GUIDE.txt` first.** This original page describes the earlier model-free viewer only.

# Additive authenticated visualization deployment

This package supplies **new files only** for an `htgaa-automation` checkout. It does not rewrite the existing GP cluster page, injector, data-fetch scripts or stylesheet.

1. Work from a fresh, authorized clone because the remote history was rewritten. Review current SvelteKit and Keycloak imports before copying any files.
2. Copy the contents of `repo_additions/` into the checkout, preserving directory structure. The new page path is `/cfps/model-free`, with private JSON at `/cfps/model-free-data`.
3. Run the notebook or CLI with an authorized source. Keep `model_free_analysis.json` in **ignored server-side** `data/cfps-model-free/`, not `static/` or Git. Cluster dataset fetching via `npm run clusters:fetch -- --token=<agent token>` is the site's *separate authorized workflow*; this package does not require or modify it.
4. Run `npm install && npm run build`. Verify: unauthorized route is gated; unauthorized `/cfps/model-free-data` returns HTTP 401; authorized GET returns schema `cfps-model-free-visualization/v1`; both PocketBase-only and published-snapshot mode menus render accurately.
5. Confirm that public page loads the published cluster stylesheet `/cfps/clusters-data/clusters.css` and the scoped companion CSS/JS without changing the original GP-informed clusters route.

**Privacy:** Recipe-level JSON contains identifying submission IDs and must remain protected. A `.gitignore` is included, but access controls and deployment secrets still require independent review against the current application checkout.

The standalone viewer can open a private JSON file via a browser file picker, without uploading its contents. It does not include that JSON in its HTML source.
