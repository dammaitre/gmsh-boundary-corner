---
name: project-state
description: BoundaryCorner implementation status, architecture decisions, open issues
metadata:
  type: project
---

## Current state (branch: complete-profile)

BoundaryCorner now covers the **full** body profile — no separate BoundaryLayer field needed.

### Parametrization change (latest)
- **Dropped**: `StartPoint` parameter and profile split
- **Added**: `NbCornerColumns` — number of arc-length-compressed columns from AxisPoint outward
- Column layout: `NbCornerColumns` compressed columns near corner + constant `MaxColumnWidth` columns for the rest of the profile
- `operator()` background field: interpolates lBLeff_→w0max_ over the compressed zone arc-length; 1e22 beyond

### Active code path: `buildForFace` (pre-mesh, in `meshGFace.cpp`)
- Orientation: uses axisPoint_ alone to orient baseVerts (back = closest to axisPoint)
- No bcStart detection needed in practice (no BL field), but kept for robustness
- N = baseVerts.size() - 2 (excludes last axis column to avoid collinear vertex issue)
- All 1D mesh nodes become column bases; normal extrusion with geometric growth

### Legacy code path: `buildCornerColumns` (post-mesh, kept for reference)
- Now finds t_start from the non-axis GVertex endpoint
- Uses nbCornerColumns_ compressed + ceil(S_rem/w0max_) constant columns
- Scales compressed section to span exactly S_total - N_constant*w0max_

### Test: test_bc.py
- Single full quarter-ellipse arc (p_top→p_nose)
- BC field only (no BL)
- Expected: ~348 quads, ~2192 triangles

**Why:** BC handles normal propagation well; profile split + BL was redundant complexity.
**How to apply:** When modifying BC, remember the full profile is in one GEdge; bcStart will be 0 without BL.
