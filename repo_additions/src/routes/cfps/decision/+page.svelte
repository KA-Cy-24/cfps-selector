<script>
  import { onMount } from 'svelte';
  import { page } from '$app/state';
  import { loginHref } from '../loginHref.js';
  let { data } = $props();
  onMount(() => {
    if (data.locked) return;
    const script = document.createElement('script');
    script.src = '/cfps/decision/decision.js';
    script.onload = () => window.CFPSDecisionLab?.mount(document.getElementById('cfps-decision-root'), {
      dataUrl: '/cfps/decision-data'
    });
    document.body.appendChild(script);
    return () => script.remove();
  });
</script>
<svelte:head>
  <title>CFPS Plate Decision Lab · HTGAA</title>
  <link rel="stylesheet" href="/cfps/decision/decision.css" />
</svelte:head>
{#if data.locked}
  <main style="background:#07101a;color:#e8eef6;padding:4rem;min-height:100vh">
    <h1>Sign in to use the CFPS Plate Decision Lab</h1>
    <p>Selection data contains private participant submissions.</p>
    <a href={loginHref(page.url)}>Sign in</a>
  </main>
{:else}
  <div id="cfps-decision-root" aria-live="polite"></div>
{/if}
