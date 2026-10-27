<script>
  // Additive companion route, does not import/modify CFPS Clusters source.
  import { onMount } from 'svelte';
  import { page } from '$app/state';
  import { loginHref } from '../loginHref.js';
  let { data } = $props();
  onMount(() => {
    if (data.locked) return;
    const script = document.createElement('script');
    script.src = '/cfps/model-free/model_free_injector.js';
    script.onload = () => {
      window.CFPSModelFreeView?.mount(document.getElementById('cfps-model-free-root'), {
        dataUrl: '/cfps/model-free-data'
      });
    };
    document.body.appendChild(script);
    return () => script.remove();
  });
</script>
<svelte:head>
  <title>CFPS Selection Coverage Analysis · HTGAA</title>
  <!-- Reuse the published CFPS cluster stylesheet and add scoped companion styles. -->
  <link rel="stylesheet" href="/cfps/clusters-data/clusters.css" />
  <link rel="stylesheet" href="/cfps/model-free/model_free.css" />
</svelte:head>
<div class="clp-page mf-page">
  <div class="clp-inner mf-inner">
    <nav class="mf-nav"><a href="/cfps">← Designer</a><a href="/cfps/clusters">GP-informed clusters</a><a href="/cfps/decision">Plate Decision Lab</a></nav>
    {#if data.locked}
      <section class="mf-locked">
        <h1>Sign in to view the CFPS model-free selection analysis</h1>
        <p>This view contains private submitted design IDs and is available only to authenticated participants.</p>
        <a href={loginHref(page.url)}>Sign in</a>
      </section>
    {:else}
      <div id="cfps-model-free-root" class="mf-root" aria-live="polite">Loading private model-free analysis…</div>
    {/if}
  </div>
</div>
