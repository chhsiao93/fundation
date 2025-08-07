# Foundation Bearing Capacity Calculation

## Overview

This document describes the methodology and equations used in the Foundation Bearing Capacity Calculator for shallow foundations in layered soil systems. The calculator implements Meyerhof's general bearing capacity equation with enhancements for layered soils using the Effective Shear Depth (He) approach.

## Theoretical Background

### Meyerhof's General Bearing Capacity Equation

The ultimate bearing capacity is calculated using Meyerhof's general equation:

```
qu = c·Nc·sc·dc·ic·ec + σ'·Nq·sq·dq·iq·eq + 0.5·γ'·B·Nγ·sγ·dγ·iγ·eγ
```

Where:
- **qu** = Ultimate bearing capacity (kPa)
- **c** = Effective cohesion (kPa)
- **σ'** = Effective overburden stress at foundation level (kPa)
- **γ'** = Effective unit weight of soil within shear zone (kN/m³)
- **B** = Foundation width (m)
- **N factors** = Bearing capacity factors (Nc, Nq, Nγ)
- **s factors** = Shape factors (sc, sq, sγ)
- **d factors** = Depth factors (dc, dq, dγ)
- **i factors** = Load inclination factors (ic, iq, iγ)
- **e factors** = Eccentricity factors (ec, eq, eγ)

### Three Components of Bearing Capacity

1. **Cohesion Term (qc)**: `c·Nc·sc·dc·ic·ec`
2. **Surcharge Term (qq)**: `σ'·Nq·sq·dq·iq·eq`
3. **Weight Term (qγ)**: `0.5·γ'·B·Nγ·sγ·dγ·iγ·eγ`

## Bearing Capacity Factors

### Basic Bearing Capacity Factors

```
Nq = e^(π·tan(φ)) · tan²(π/4 + φ/2)

Nc = (Nq - 1) / tan(φ)    for φ > 0
Nc = 5.7                  for φ = 0

Nγ = (Nq - 1) · tan(1.4φ)  (Meyerhof's formula)
```

Where φ is the internal friction angle in radians.

### Shape Factors

For rectangular foundations with length L and width B:

```
B/L = width-to-length ratio
Kp = tan²(π/4 + φ/2)  (passive earth pressure coefficient)

sc = 1 + 0.2·Kp·(B/L)
sq = 1 + 0.1·Kp·(B/L)    for φ > 10°
sq = 1.0                 for φ ≤ 10°
sγ = sq
```

### Depth Factors

For foundations with depth D and width B:

```
D/B = depth-to-width ratio
Kp = tan²(π/4 + φ/2)

dc = 1 + 0.2·√(Kp)·(D/B)
dq = 1 + 0.1·√(Kp)·(D/B)    for φ > 10°
dq = 1.0                     for φ ≤ 10°
dγ = dq
```

### Load Inclination Factors

For inclined loads with angle λ (degrees):

```
ic = (1 - λ/90)²
iq = (1 - λ/90)²
iγ = (1 - λ/φ)²    for φ > 10°
iγ = 1.0           for φ ≤ 10°
```

**Note**: For inclined loads (λ > 0), shape factors are set to 1.0.

### Eccentricity Factors

For eccentric loads with eccentricity e:

```
ec = 1 - 2·(e/B)
eq = 1 - 2·(e/B)
eγ = 1 - 3.5·(e/B) + 3·(e/B)²
```

Where e < B/2 (eccentricity must be less than half the foundation width).

## Layered Soil Analysis

### Effective Shear Depth (He) Method

For layered soils, the calculator uses the Effective Shear Depth approach with iterative convergence:

#### Initial Calculation
```
He₀ = 0.5 · B · tan(45° + φ₀/2)
```

Where φ₀ is the friction angle of the soil layer immediately below the foundation.

#### Iterative Process

1. **Calculate weighted average properties** within current He:
   ```
   cav = Σ[ci · ti] / Σ[ti]
   φav = arctan(Σ[tan(φi) · ti] / Σ[ti])
   ```
   
   Where ci, φi, ti are cohesion, friction angle, and thickness contribution of layer i within He.

2. **Recalculate He** using average friction angle:
   ```
   He(new) = 0.5 · B · tan(45° + φav/2)
   ```

3. **Check convergence**:
   ```
   |He(new) - He(old)| < tolerance (default: 0.0001 m)
   ```

4. **Repeat** until convergence (typically 1-10 iterations).

### Overburden Effective Stress Calculation

For layered soils with groundwater:

```
σ' = Σ[γi · hi] - uw

where:
uw = γw · (depth - GWT depth)  if depth > GWT depth
uw = 0                         if depth ≤ GWT depth
```

Where:
- γi = unit weight of layer i
- hi = thickness of layer i contributing to the depth
- γw = unit weight of water (9.81 kN/m³)

### Average Effective Unit Weight for qγ Term

Using the downward-tapering wedge formula:

```
γe = (2/He²) · Σ[γi' · ((zi+1-zi)·He - 0.5·(zi+1² - zi²))]

where:
γi' = γi - γw  (effective unit weight if below GWT)
γi' = γi       (total unit weight if above GWT)
zi, zi+1 = layer boundaries from foundation level
```

## Calculation Workflow

### 1. Input Processing
- Load configuration from YAML file
- Validate input parameters
- Check foundation type and dimensions

### 2. Layered Soil Analysis
- Calculate Effective Shear Depth (He) with iterative convergence
- Determine weighted average cohesion (cav) and friction angle (φav)
- Calculate overburden effective stress at foundation level

### 3. Factor Calculations
- Compute bearing capacity factors (Nc, Nq, Nγ) using φav
- Calculate shape factors based on foundation geometry
- Determine depth factors based on D/B ratio
- Apply load inclination factors if λ > 0
- Apply eccentricity factors if e > 0

### 4. Unit Weight Analysis
- Calculate average effective unit weight within He using wedge formula
- Account for groundwater effects (buoyancy)

### 5. Ultimate Bearing Capacity
```
qu = cav·Nc·sc·dc·ic·ec + σ'·Nq·sq·dq·iq·eq + 0.5·γe·B·Nγ·sγ·dγ·iγ·eγ
```

### 6. Allowable Bearing Capacity
```
qa = qu / FS
```

Where FS is the factor of safety (typically 2.5-3.0).

## Key Features

### Layered Soil Handling
- Supports unlimited number of soil layers
- Automatic layer boundary detection within shear zone
- Weighted averaging based on layer contributions

### Groundwater Effects
- Effective stress calculations with pore pressure
- Buoyancy effects on unit weight below GWT
- Automatic handling of partially submerged layers

### Foundation Types
- Strip, square, circular, and rectangular foundations
- Automatic shape factor selection
- Dimension validation for rectangular foundations

### Loading Conditions
- Vertical and inclined loads
- Eccentric loading with automatic factor adjustments
- Surcharge loads at ground surface

### Advanced Features
- Iterative convergence for He calculation
- Professional soil profile diagram generation
- Comprehensive result breakdown by components
- Input validation and error handling

## Example Calculation

For a 2m × 2m square foundation at 1.5m depth with:
- Applied load: 500 kPa
- Layered soil with different φ values
- Groundwater at 2.0m depth

**Results:**
- He = 1.73 m (converged in 1 iteration)
- φav = 30.0°
- Ultimate bearing capacity = 2245 kPa
- Allowable capacity (FS=3) = 748 kPa

## References

1. Meyerhof, G.G. (1963). "Some recent research on the bearing capacity of foundations." Canadian Geotechnical Journal, 1(1), 16-26.

2. Das, B.M. (2019). "Principles of Foundation Engineering." 9th Edition, Cengage Learning.

3. Bowles, J.E. (1997). "Foundation Analysis and Design." 5th Edition, McGraw-Hill.

4. Coduto, D.P., Yeung, M.R., and Kitch, W.A. (2016). "Geotechnical Engineering: Principles and Practices." 2nd Edition, Pearson.