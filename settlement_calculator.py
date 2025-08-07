#!/usr/bin/env python3
"""
Foundation Settlement Calculator

This script calculates the settlement of shallow foundations using consolidation theory
for layered soil systems with the same configuration format as the bearing capacity calculator.

Author: Foundation Analysis Tool
"""

import math
import yaml
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional
from enum import Enum


class SettlementConfig:
    """Configuration class for settlement analysis using the same format as bearing capacity."""
    
    def __init__(self, config_path: str = "config.yaml"):
        """
        Initialize configuration from YAML file (same as bearing capacity config).
        
        Args:
            config_path: Path to the YAML configuration file
        """
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self._validate_config()
    
    def _load_config(self) -> Dict:
        """Load configuration from YAML file."""
        if not self.config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {self.config_path}")
        
        try:
            with open(self.config_path, 'r') as file:
                config = yaml.safe_load(file)
                return config
        except yaml.YAMLError as e:
            raise ValueError(f"Error parsing YAML configuration: {e}")
    
    def _validate_config(self) -> None:
        """Validate the loaded configuration."""
        required_sections = ['soil', 'foundation']
        for section in required_sections:
            if section not in self.config:
                raise ValueError(f"Missing required section in config: {section}")
    
    # Soil Properties (same format as bearing capacity)
    @property
    def soil_cohesion(self) -> List[float]:
        return self.config['soil']['cohesion']
    
    @property
    def soil_unit_weight(self) -> List[float]:
        return self.config['soil']['unit_weight']
    
    @property
    def soil_friction_angle(self) -> List[float]:
        return self.config['soil']['friction_angle']
    
    @property
    def soil_thickness(self) -> List[float]:
        return self.config['soil']['thickness']
    
    @property
    def compression_index(self) -> List[float]:
        return self.config['soil']['compression_index']
    
    @property
    def recompression_index(self) -> List[float]:
        return self.config['soil']['recompression_index']
    
    @property
    def void_ratio(self) -> List[float]:
        return self.config['soil']['void_ratio']

    @property
    def ocr(self) -> List[float]:
        return self.config['soil']['ocr']

    @property
    def gwt_depth(self) -> Optional[float]:
        return self.config['soil'].get('gwt')
    
    @property
    def foundation_width(self) -> float:
        return self.config['foundation']['width']
    
    @property
    def foundation_length(self) -> Optional[float]:
        return self.config['foundation']['length']
    
    @property
    def foundation_depth(self) -> float:
        return self.config['foundation']['depth']
    
    # Loading (from config or parameter)
    @property
    def applied_load(self) -> float:
        return self.config['loading']['load']
    
    # Output Settings
    @property
    def show_components(self) -> bool:
        return self.config['output'].get('show_components', True)
    
    @property
    def decimal_places(self) -> int:
        return self.config['output'].get('decimal_places', 3)


class SettlementCalculator:
    def __init__(self):
        self.gamma_w = 9.81  # Water unit weight (kN/m³)
    
    def _validate_positive_number(self, value: float, name: str) -> None:
        """Validate that a value is a positive number."""
        if not isinstance(value, (int, float)) or value <= 0:
            raise ValueError(f"{name} must be a positive number, got {value}")
    
    def _calculate_effective_stress(self, unit_weights: List[float], 
                                  layer_thicknesses: List[float],
                                  gwt_depth: Optional[float],
                                  target_depth: float) -> float:
        """
        Calculate effective stress at a given depth.
        
        Args:
            unit_weights: List of unit weights for each layer (kN/m³)
            layer_thicknesses: List of thicknesses for each layer (m)
            gwt_depth: Groundwater table depth from surface (m)
            target_depth: Target depth from surface (m)
            
        Returns:
            Effective stress at target depth (kPa)
        """
        total_stress = 0.0
        pore_pressure = 0.0
        current_depth = 0.0
        
        # Calculate total vertical stress
        for gamma, thickness in zip(unit_weights, layer_thicknesses):
            layer_bottom = current_depth + thickness
            
            if current_depth >= target_depth:
                break
            
            # Determine effective thickness within target depth
            effective_thickness = min(thickness, target_depth - current_depth)
            
            # Add to total stress
            total_stress += gamma * effective_thickness
            current_depth += thickness
            
        # Calculate pore pressure if below GWT
        pore_pressure = 0.0
        if gwt_depth is not None and target_depth > gwt_depth:
            water_depth = target_depth - gwt_depth
            pore_pressure = self.gamma_w * water_depth
            
        
        # Effective stress = Total stress - Pore pressure
        return max(0.0, total_stress - pore_pressure)
    
    def calculate_consolidation_settlement(self,
                                         applied_load: float,
                                         unit_weights: List[float],
                                         layer_thicknesses: List[float],
                                         compression_indices: List[float],
                                         recompression_indices: List[float],
                                         void_ratios: List[float],
                                         ocr: List[float],
                                         gwt_depth: Optional[float],
                                         foundation_width: float,
                                         foundation_length: float,
                                         foundation_depth: float) -> Dict[str, float]:
        """
        Calculate consolidation settlement for layered soil using same config format as bearing capacity.
        Args:
            applied_load: Applied pressure (kPa)
            unit_weights: List of unit weights for each layer (kN/m³)
            layer_thicknesses: List of thicknesses for each layer (m)
            compression_indices: List of compression indices for each layer
            void_ratios: List of initial void ratios for each layer
            preconsolidation_pressures: List of preconsolidation pressures for each layer (kPa)
            gwt_depth: Groundwater table depth from surface (m)
            foundation_width: Foundation width (m)
            foundation_length: Foundation length (m)
            foundation_depth: Foundation depth (m)
            
        Returns:
            Dictionary with settlement components
        """
        applied_pressure = applied_load / foundation_width / foundation_length
        total_settlement = 0.0
        layer_settlements = []
        
        current_depth_from_surface = 0.0
        analysis_bottom = sum(layer_thicknesses)
        
        for i, (thickness, gamma, cc, cr, e0, ocr) in enumerate(zip(
            layer_thicknesses, unit_weights, compression_indices, recompression_indices, void_ratios, ocr)):
            layer_top = current_depth_from_surface
            layer_bottom = current_depth_from_surface + thickness
            
            # Skip if layer is completely above foundation
            if layer_bottom <= foundation_depth:
                current_depth_from_surface += thickness
                continue
            
            # Skip if layer starts below analysis depth
            if layer_top >= analysis_bottom:
                break
            
            # Adjust layer boundaries to analysis zone
            effective_layer_top = max(layer_top, foundation_depth)
            effective_layer_bottom = min(layer_bottom, analysis_bottom)
            effective_thickness = effective_layer_bottom - effective_layer_top
            
            if effective_thickness <= 0:
                current_depth_from_surface += thickness
                continue
            
            # Calculate stress increase at mid-depth of effective layer
            mid_depth_from_foundation = (effective_layer_top + effective_layer_bottom) / 2 - foundation_depth
            stress_increase = self._calculate_stress_increase(
                applied_pressure, foundation_width, mid_depth_from_foundation
            )
            
            # Calculate initial effective stress at mid-depth
            mid_depth_from_surface = foundation_depth + mid_depth_from_foundation
            initial_effective_stress = self._calculate_effective_stress(
                unit_weights, layer_thicknesses, gwt_depth, mid_depth_from_surface
            )
            pc = initial_effective_stress * ocr  # Preconsolidation pressure
            # Calculate consolidation settlement using Terzaghi's theory
            # S = (Cc * H * log((σ'₀ + Δσ) / σ'₀)) / (1 + e₀) for normally consolidated soil
            # S = (Cc * H * log((σ'₀ + Δσ) / σ'pc)) / (1 + e₀) for overconsolidated soil where σ'₀ + Δσ > σ'pc
            if initial_effective_stress + stress_increase <= pc:
                # Overconsolidated - use recompression index (typically Cc/5 to Cc/10, use Cc/7)
                layer_settlement = (cr * effective_thickness * 
                                  math.log10((initial_effective_stress + stress_increase) / 
                                           initial_effective_stress)) / (1 + e0)
            else:
                if initial_effective_stress <= pc:
                    # Transition from overconsolidated to normally consolidated
                    settlement_oc = (cr * effective_thickness * 
                                   math.log10(pc / initial_effective_stress)) / (1 + e0)
                    settlement_nc = (cc * effective_thickness * 
                                   math.log10((initial_effective_stress + stress_increase) / pc)) / (1 + e0)
                    layer_settlement = settlement_oc + settlement_nc
                else:
                    # Normally consolidated
                    layer_settlement = (cc * effective_thickness * 
                                      math.log10((initial_effective_stress + stress_increase) / 
                                               initial_effective_stress)) / (1 + e0)
            
            total_settlement += layer_settlement
            layer_settlements.append({
                'layer_index': i + 1,
                'depth_range': f"{effective_layer_top:.1f}-{effective_layer_bottom:.1f}m",
                'thickness': effective_thickness,
                'initial_stress': initial_effective_stress,
                'stress_increase': stress_increase,
                'final_stress': initial_effective_stress + stress_increase,
                'precon_pressure': pc,
                'settlement': layer_settlement
            })
            
            current_depth_from_surface += thickness
        
        return {
            'total_settlement': total_settlement,
            'layer_settlements': layer_settlements,
            'method': 'Consolidation Theory (Terzaghi)'
        }
    
    def _calculate_stress_increase(self, applied_pressure: float, 
                                 foundation_width: float, depth: float) -> float:
        """
        Calculate stress increase using 2:1 method (simplified).
        More sophisticated methods like Boussinesq can be implemented.
        """
        if depth <= 0:
            return applied_pressure
        
        # 2:1 load distribution method
        # Stress spreads outward at 2 horizontal : 1 vertical ratio
        loaded_area_at_depth = (foundation_width + depth) ** 2
        foundation_area = foundation_width ** 2
        
        stress_reduction_factor = foundation_area / loaded_area_at_depth
        
        return applied_pressure * stress_reduction_factor


def run_settlement_analysis(config_path: str) -> None:
    """
    Run settlement analysis using configuration file.
    
    Args:
        config_path: Path to the YAML configuration file
    """
    try:
        # Load configuration
        config = SettlementConfig(config_path)
        calculator = SettlementCalculator()
        
        print("Foundation Settlement Calculator")
        print("=" * 60)
        print(f"Configuration file: {config.config_path}")
        print(f"Foundation: {config.foundation_width:.1f}m x {config.foundation_length or config.foundation_width:.1f}m")
        print(f"Foundation Depth: {config.foundation_depth:.1f}m")
        if config.gwt_depth is not None:
            print(f"Groundwater Table: {config.gwt_depth:.1f}m")
        print("=" * 60)
        
        # Calculate consolidation settlement
        result = calculator.calculate_consolidation_settlement(
            applied_load=config.applied_load,
            unit_weights=config.soil_unit_weight,
            layer_thicknesses=config.soil_thickness,
            compression_indices=config.compression_index,
            recompression_indices=config.recompression_index,
            void_ratios=config.void_ratio,
            ocr=config.ocr,
            gwt_depth=config.gwt_depth,
            foundation_width=config.foundation_width,
            foundation_length=config.foundation_length,
            foundation_depth=config.foundation_depth,
        )
        
        # Display results
        print(f"\n{result['method']}:")
        print("-" * 40)
        settlement = result['total_settlement']
        print(f"Total Consolidation Settlement: {settlement*1000:.1f} mm ({settlement:.4f} m)")
        
        if config.show_components and result['layer_settlements']:
            print("\nLayer-by-layer breakdown:")
            print(f"{'Layer':^6} | {'Depth Range':^12} | {'σ₀ (kPa)':^8} | {'Δσ (kPa)':^8} | {'σpc (kPa)':^9} | {'Settlement':^10}")
            print("-" * 70)
            for layer_result in result['layer_settlements']:
                print(f"{layer_result['layer_index']:^6} | "
                      f"{layer_result['depth_range']:^12} | "
                      f"{layer_result['initial_stress']:^8.1f} | "
                      f"{layer_result['stress_increase']:^8.1f} | "
                      f"{layer_result['precon_pressure']:^9.1f} | "
                      f"{layer_result['settlement']*1000:^8.1f} mm")
        
        print("\n" + "=" * 60)
        print("SUMMARY:")
        print("=" * 60)
        print(f"Total Settlement: {settlement*1000:.1f} mm ({settlement:.4f} m)")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Foundation Settlement Calculator (Consolidation Method Only)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python settlement_calculator.py                    # Calculate settlement
  python settlement_calculator.py -c config.yaml     # Use custom config file
        """
    )
    
    parser.add_argument(
        '-c', '--config',
        default='config.yaml',
        help='Path to YAML configuration file (default: config.yaml - same as bearing capacity)'
    )
    
    args = parser.parse_args()
    
    # Run settlement analysis with configuration
    run_settlement_analysis(args.config)


if __name__ == "__main__":
    main()