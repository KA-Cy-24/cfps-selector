// Gated JSON: NEVER put model_free_analysis.json under static/ or commit it.
// Separate endpoint to protect private recipe IDs on direct URL access.
import { error } from '@sveltejs/kit';
import { requestUser } from '$lib/server/auth.js';
import { readFile } from 'node:fs/promises';
import { resolve } from 'node:path';
export async function GET(event) {
  const user = await requestUser(event);
  if (!user) throw error(401, 'Sign in to access CFPS selection analysis.');
  const file = resolve(process.cwd(), 'data', 'cfps-model-free', 'model_free_analysis.json');
  let raw;
  try { raw = await readFile(file); }
  catch (e) {
    if (e?.code === 'ENOENT') throw error(503, 'Model-free analysis has not been generated on this server.');
    throw error(503, 'Private analysis temporarily unavailable.');
  }
  return new Response(raw, {
    headers: {
      'Content-Type': 'application/json; charset=utf-8',
      'Cache-Control': 'private, no-store',
      'X-Content-Type-Options': 'nosniff',
      'Content-Security-Policy': "default-src 'none'",
    }
  });
}
