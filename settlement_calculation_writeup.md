# Foundation Settlement Calculation

## Overview

This document describes the methodology and equations used in the Foundation Settlement Calculator for shallow foundations in layered soil systems. The calculator implements consolidation settlement analysis using Terzaghi's consolidation theory with enhancements for layered soils, overconsolidated behavior, and groundwater effects.

## Theoretical Background

### Consolidation Settlement Theory

Settlement in cohesive soils occurs primarily due to consolidation - the gradual expulsion of pore water under applied loads. The settlement calculation is based on Terzaghi's one-dimensional consolidation theory.

### Total Settlement Components

Total settlement typically consists of three components:

1. **Immediate (Elastic) Settlement**: Occurs instantly upon load application
2. **Primary Consolidation Settlement**: Time-dependent settlement due to pore water expulsion
3. **Secondary Compression**: Long-term settlement under constant effective stress

*Note: This calculator focuses on primary consolidation settlement only.*

## Primary Consolidation Settlement

### Basic Consolidation Equation

The primary consolidation settlement for a soil layer is given by:

```
S = (Cc × H × log₁₀((σ'₀ + Δσ) / σ'₀)) / (1 + e₀)
```

Where:
- **S** = Settlement of the layer (m)
- **Cc** = Compression index (dimensionless)
- **H** = Layer thickness (m)
- **σ'₀** = Initial effective stress at mid-depth of layer (kPa)
- **Δσ** = Stress increase due to applied load (kPa)
- **e₀** = Initial void ratio (dimensionless)

### Overconsolidated Soil Behavior

For overconsolidated soils, the settlement calculation depends on the relationship between current stress and preconsolidation pressure:

#### Case 1: Stress remains below preconsolidation pressure
```
σ'₀ + Δσ ≤ σ'pc

S = (Cr × H × log₁₀((σ'₀ + Δσ) / σ'₀)) / (1 + e₀)
```

#### Case 2: Stress exceeds preconsolidation pressure
```
σ'₀ + Δσ > σ'pc

S = S₁ + S₂

where:
S₁ = (Cr × H × log₁₀(σ'pc / σ'₀)) / (1 + e₀)        [Recompression]
S₂ = (Cc × H × log₁₀((σ'₀ + Δσ) / σ'pc)) / (1 + e₀)  [Virgin compression]
```

Where:
- **σ'pc** = Preconsolidation pressure (kPa)
- **Cr** = Recompression index (typically Cc/5 to Cc/10)
- **OCR** = Overconsolidation ratio = σ'pc / σ'₀

### Preconsolidation Pressure

The preconsolidation pressure is calculated from the overconsolidation ratio:

```
σ'pc = OCR × σ'₀
```

## Stress Calculations

### Initial Effective Stress

For layered soil systems with groundwater, the initial effective stress at any depth is:

```
σ'₀ = Σ[γᵢ × hᵢ] - uw

where:
uw = γw × (depth - GWT depth)  if depth > GWT depth
uw = 0                         if depth ≤ GWT depth
```

Where:
- **γᵢ** = Unit weight of layer i (kN/m³)
- **hᵢ** = Thickness of layer i contributing to the depth (m)
- **γw** = Unit weight of water (9.81 kN/m³)
- **uw** = Pore water pressure (kPa)

### Stress Increase Due to Applied Load

The stress increase at depth is calculated using the 2:1 load distribution method:

```
Δσ = q × (Bf²) / (Bf + z)²

where:
q = Applied pressure = P / (B × L)
Bf = Foundation width (m)
z = Depth below foundation level (m)
```

*Note: More sophisticated methods like Boussinesq or Westergaard equations can be implemented for improved accuracy.*

## Layered Soil Analysis

### Layer-by-Layer Settlement Calculation

For layered soil systems, settlement is calculated for each compressible layer and summed:

```
Stotal = Σ Si
```

### Layer Contribution Analysis

Each layer's contribution is evaluated based on:

1. **Layer boundaries**: Only the portion below foundation level is considered
2. **Effective thickness**: Intersection with the analysis depth zone
3. **Mid-depth properties**: Stress calculations at layer mid-depth
4. **Consolidation parameters**: Layer-specific Cc, Cr, e₀, OCR values

### Analysis Depth

The analysis typically extends to a depth where stress increase becomes negligible:

```
Analysis depth = Foundation depth + Σ(layer thicknesses)
```

Or until stress increase < 10% of initial effective stress.

## Calculation Workflow

### 1. Input Processing
- Load configuration from YAML file (same format as bearing capacity)
- Validate soil parameters and foundation geometry
- Process applied load and foundation dimensions

### 2. Foundation Pressure Calculation
```
Applied pressure (q) = Applied load (P) / Foundation area (B × L)
```

### 3. Layer Analysis Loop
For each soil layer below foundation:

#### a. Determine Effective Layer Boundaries
```
Effective top = max(layer top, foundation depth)
Effective bottom = min(layer bottom, analysis depth)
Effective thickness = Effective bottom - Effective top
```

#### b. Calculate Mid-Depth Properties
```
Mid-depth from surface = (Effective top + Effective bottom) / 2
Mid-depth from foundation = Mid-depth from surface - Foundation depth
```

#### c. Initial Effective Stress
Calculate σ'₀ at mid-depth using layered soil effective stress formula.

#### d. Stress Increase
Calculate Δσ at mid-depth using 2:1 method or other distribution.

#### e. Preconsolidation Pressure
```
σ'pc = OCR × σ'₀
```

#### f. Settlement Calculation
Apply appropriate consolidation equation based on stress state relative to preconsolidation pressure.

### 4. Total Settlement
Sum all layer settlements:
```
Total Settlement = Σ(Layer settlements)
```

## Key Parameters

### Soil Properties (from config.yaml)
```yaml
soil:
  compression_index: [Cc1, Cc2, ...]     # Compression index for each layer
  recompression_index: [Cr1, Cr2, ...]   # Recompression index for each layer  
  void_ratio: [e01, e02, ...]            # Initial void ratio for each layer
  ocr: [OCR1, OCR2, ...]                 # Overconsolidation ratio for each layer
  unit_weight: [γ1, γ2, ...]             # Unit weight for each layer (kN/m³)
  thickness: [H1, H2, ...]               # Layer thickness (m)
  gwt: 2.0                               # Groundwater table depth (m)
```

### Typical Parameter Values

| Parameter | Symbol | Typical Range | Units |
|-----------|---------|---------------|-------|
| Compression Index | Cc | 0.1 - 1.5 | - |
| Recompression Index | Cr | Cc/5 - Cc/10 | - |
| Initial Void Ratio | e₀ | 0.4 - 2.0 | - |
| Overconsolidation Ratio | OCR | 1.0 - 10+ | - |
| Unit Weight | γ | 16 - 22 | kN/m³ |

## Settlement Analysis Features

### Layered Soil Handling
- Supports unlimited number of soil layers
- Individual consolidation parameters for each layer
- Automatic layer boundary detection within analysis zone
- Proper treatment of partial layers

### Groundwater Effects
- Effective stress calculations with pore pressure
- Buoyancy effects on unit weight below groundwater table
- Automatic handling of partially submerged layers

### Overconsolidation Analysis
- OCR-based preconsolidation pressure calculation
- Recompression vs. virgin compression behavior
- Transition handling for mixed stress states

### Loading Conditions
- Applied load converted to foundation pressure
- 2:1 stress distribution method (expandable to other methods)
- Depth-dependent stress reduction

### Output Features
- Total settlement in mm and meters
- Layer-by-layer breakdown with stress states
- Initial stress, stress increase, and final stress for each layer
- Preconsolidation pressure comparison
- Settlement contribution per layer

## Example Calculation

For a 2m × 2m square foundation with applied load of 1000 kN:

**Input Parameters:**
- Applied pressure: 1000 kN / (2m × 2m) = 250 kPa
- Foundation depth: 1.0 m
- Layered soil with different consolidation parameters
- Groundwater at 2.0 m depth

**Layer 1 (1.0-3.0m depth):**
- Thickness: 2.0 m (effective thickness: 2.0 m)
- Mid-depth: 2.0 m from surface, 1.0 m from foundation
- Initial effective stress: 18.0 kPa
- Stress increase: 125 kPa (2:1 method)
- OCR: 2.0 → σ'pc = 36.0 kPa
- Since σ'₀ + Δσ = 143 kPa > σ'pc: Mixed behavior
- Settlement: 15.2 mm

**Layer 2 (3.0-18.0m depth):**
- Effective thickness: 15.0 m
- Settlement contribution: 42.8 mm

**Results:**
- Total consolidation settlement: 58.0 mm
- Allowable settlement typically: 25-50 mm for buildings

## Limitations and Considerations

### Method Limitations
1. **One-dimensional consolidation**: Assumes vertical drainage only
2. **Linear stress-strain relationship**: Uses semi-logarithmic compression
3. **Simplified stress distribution**: 2:1 method is approximate
4. **No time-rate analysis**: Provides ultimate settlement only

### Design Considerations
1. **Allowable settlement**: Typically 25-50 mm for buildings
2. **Differential settlement**: Often more critical than total settlement
3. **Construction sequence**: Staged loading effects not considered
4. **Secondary compression**: Long-term settlement beyond primary consolidation

### Accuracy Improvements
1. **Boussinesq stress distribution**: More accurate than 2:1 method
2. **Non-linear consolidation**: For heavily overconsolidated soils
3. **Multi-dimensional analysis**: For complex loading and geometry
4. **Time-rate analysis**: Using Terzaghi's time-consolidation theory

## References

1. Terzaghi, K. (1943). "Theoretical Soil Mechanics." John Wiley & Sons, New York.

2. Das, B.M. (2019). "Principles of Foundation Engineering." 9th Edition, Cengage Learning.

3. Holtz, R.D., Kovacs, W.D., and Sheahan, T.C. (2011). "An Introduction to Geotechnical Engineering." 2nd Edition, Pearson.

4. Coduto, D.P., Yeung, M.R., and Kitch, W.A. (2016). "Geotechnical Engineering: Principles and Practices." 2nd Edition, Pearson.

5. Bowles, J.E. (1997). "Foundation Analysis and Design." 5th Edition, McGraw-Hill.

6. Duncan, J.M. and Buchignani, A.L. (1976). "An Engineering Manual for Settlement Studies." Department of Civil Engineering, University of California, Berkeley.