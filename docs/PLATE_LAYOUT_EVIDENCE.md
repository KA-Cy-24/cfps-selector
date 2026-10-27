# Plate geometry and sample replication evidence

Inspected source: authorized supplied `cfps-plates.zip` archive (external plate designs/results, not eligible submitted candidates). The examined Ginkgo-style plate records specify:

| Property | Inspected value | Consequence |
|---|---|---|
| `n_rows` | 16 | Rows A–P |
| `n_columns` | 24 | Columns 1–24 |
| `reserved_columns` | `[1,2,3,4]` | 64 reserved wells, unavailable for participant-design allocation |
| Experimental columns | 5–24 | 320 sample positions |
| `replicate_factor` | 4 | Normally 80 distinct experimental conditions, each repeated 4 times |
| `array_order` | `column` | The illustrative well map fills column-major; laboratory may randomize |

Example plates used to cross-check: `styczynski_sfgfp_p1.json`, `styczynski_petase_p3.json`, and `adamala_lab_PETase_Design_v3_validated.json`. In examined sample arrays, all 320 samples are marked `experimental`; 80 distinct sample IDs occur exactly four times. Other supplied example plates also publish the same reserved-columns/replication metadata.

Externally verified source: [Ginkgo CFPS automation repository](https://github.com/ginkgobioworks/ginkgo-automation-cfps) and [plate-model explanation](https://deepwiki.com/ginkgobioworks/ginkgo-automation-cfps/4.3-plate-model). This is the **Ginkgo automation CFPS plate schema**, **not** the different commercial CFPS Ready384 reagent-kit layout. These are not interchangeable.

**Warning:** The final competition execution specification has not been independently verified. Do not treat the displayed mapping as robot executable, and have the operator confirm replicate counts and control placement. The UI makes technical replication explicit so a 320-submission plan cannot quietly be misrepresented as a 320-unique fourfold-replicate physical plate.
