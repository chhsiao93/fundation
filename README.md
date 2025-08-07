# Foundation Analysis Tools

A comprehensive Python toolkit for foundation engineering analysis, including bearing capacity and settlement calculations for shallow foundations in layered soil systems.

## Features

### 🏗️ Bearing Capacity Analysis
- **Meyerhof's General Equation**: Complete implementation with all bearing capacity factors
- **Layered Soil Support**: Effective Shear Depth (He) method with iterative convergence
- **Foundation Types**: Strip, square, circular, and rectangular foundations
- **Loading Conditions**: Vertical, inclined, and eccentric loads
- **Groundwater Effects**: Automatic effective stress and buoyancy calculations

### 📐 Settlement Analysis
- **Consolidation Theory**: Terzaghi's one-dimensional consolidation settlement
- **Overconsolidated Soils**: OCR-based preconsolidation pressure handling
- **Layered Analysis**: Layer-by-layer settlement computation
- **Stress Distribution**: 2:1 method with expandable architecture

### 📊 Visual Tools
- **Soil Profile Diagrams**: ASCII-based foundation and soil layer visualization
- **Detailed Output**: Component breakdown and layer-by-layer analysis
- **Professional Reports**: Comprehensive calculation summaries

## Installation

### Prerequisites
- Python 3.7 or higher
- PyYAML library

### Setup
```bash
git clone <repository-url>
cd fundation
pip install -r requirements.txt
```

## Quick Start

### 1. Configure Your Project
Edit `config.yaml` with your foundation and soil parameters:

```yaml
soil:
  cohesion: [20.0, 25.0]          # kPa for each layer
  friction_angle: [30.0, 32.0]    # degrees for each layer
  unit_weight: [18.0, 20.0]       # kN/m³ for each layer
  thickness: [2.0, 15.0]          # m for each layer
  gwt: 2.0                        # Groundwater depth (m)
  
  # Settlement parameters
  compression_index: [0.3, 0.4]
  recompression_index: [0.05, 0.05]
  void_ratio: [0.6, 0.8]
  ocr: [2.0, 1.0]

foundation:
  width: 2.0                      # m
  length: 2.0                     # m
  depth: 1.0                      # m

loading:
  load: 1000.0                    # kN
  load_angle: 0.0                 # degrees
  eccentricity: 0.0               # m
```

### 2. Run Bearing Capacity Analysis
```bash
python bearing_capacity.py
```

### 3. Run Settlement Analysis
```bash
python settlement_calculator.py
```

## Usage Examples

### Bearing Capacity Calculation
```python
from bearing_capacity import FoundationConfig, BearingCapacityCalculator

# Load configuration
config = FoundationConfig("config.yaml")
calculator = BearingCapacityCalculator()

# Calculate bearing capacity
result = calculator.calculate_bearing_capacity(config)

print(f"Ultimate Bearing Capacity: {result['ultimate_capacity']:.1f} kPa")
print(f"Allowable Bearing Capacity: {result['allowable_capacity']:.1f} kPa")
```

### Settlement Calculation
```python
from settlement_calculator import SettlementConfig, SettlementCalculator

# Load configuration
config = SettlementConfig("config.yaml")
calculator = SettlementCalculator()

# Calculate settlement
result = calculator.calculate_consolidation_settlement(
    applied_load=config.applied_load,
    unit_weights=config.soil_unit_weight,
    layer_thicknesses=config.soil_thickness,
    # ... other parameters
)

print(f"Total Settlement: {result['total_settlement']*1000:.1f} mm")
```

## Key Algorithms

### Bearing Capacity
- **Meyerhof's Equation**: `qu = c·Nc·sc·dc·ic·ec + σ'·Nq·sq·dq·iq·eq + 0.5·γ'·B·Nγ·sγ·dγ·iγ·eγ`
- **Effective Shear Depth**: Iterative calculation using `He = 0.5·B·tan(45° + φ/2)`
- **Layered Analysis**: Weighted averaging within shear zone

### Settlement
- **Consolidation Settlement**: `S = (Cc·H·log₁₀((σ'₀ + Δσ) / σ'₀)) / (1 + e₀)`
- **Overconsolidation**: OCR-based preconsolidation pressure handling
- **Stress Distribution**: 2:1 load distribution method

## Output Examples

### Bearing Capacity Results
```
Foundation Bearing Capacity Calculator
======================================
Ultimate Bearing Capacity: 2245.3 kPa
Allowable Bearing Capacity: 748.4 kPa
Safety Factor: 3.0

Components:
- Cohesion term (qc): 1045.2 kPa
- Surcharge term (qq): 652.8 kPa  
- Weight term (qγ): 547.3 kPa
```

### Settlement Analysis Results
```
Foundation Settlement Calculator
================================
Total Consolidation Settlement: 45.2 mm (0.0452 m)

Layer-by-layer breakdown:
Layer | Depth Range | σ₀ (kPa) | Δσ (kPa) | σpc (kPa) | Settlement
  1   |  1.0-3.0m   |   18.0   |   125.0  |   36.0    |  15.2 mm
  2   |  3.0-18.0m  |   56.2   |    78.5  |   56.2    |  30.0 mm
```

## Technical Documentation

Comprehensive technical documentation is available:

- **[Bearing Capacity Analysis](bearing_capacity_writeup.md)**: Complete methodology, equations, and workflow
- **[Settlement Calculation](settlement_calculation_writeup.md)**: Consolidation theory, layered analysis, and implementation
