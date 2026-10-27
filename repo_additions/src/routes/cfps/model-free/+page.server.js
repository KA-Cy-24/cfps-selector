// Additive authenticated selection-analysis route; GP cluster source files unchanged.
import { readSession, KC_COOKIE } from '$lib/server/session.js';
import { isConfigured as kcIsConfigured } from '$lib/server/keycloak.js';
export function load({ cookies }) {
  const user = readSession(cookies.get(KC_COOKIE));
  if (!user) return { locked: true, kcConfigured: kcIsConfigured() };
  return { locked: false, kcUser: { name: user.name } };
}
