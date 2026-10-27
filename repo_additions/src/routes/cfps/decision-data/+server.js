// Protected participant-level inputs; the JSON MUST live outside tracked/static folders.
import { error } from '@sveltejs/kit';
import { readSession, KC_COOKIE } from '$lib/server/session.js';
import { readFile } from 'node:fs/promises';
import { resolve } from 'node:path';
export async function GET({ cookies }) {
  if (!readSession(cookies.get(KC_COOKIE))) throw error(401, 'Authentication required.');
  const file = resolve(process.cwd(), 'data', 'cfps-decision', 'cfps_decision_input_private.json');
  let raw;
  try { raw = await readFile(file, 'utf8'); }
  catch (e) { throw error(503, 'Decision input unavailable; generate a private input using the Colab notebook.'); }
  return new Response(raw, { headers: {
    'Content-Type': 'application/json; charset=utf-8',
    'Cache-Control': 'private, no-store',
    'X-Content-Type-Options': 'nosniff',
    'Content-Security-Policy': "default-src 'none'",
  }});
}
