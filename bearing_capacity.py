#!/usr/bin/env python3
"""
Foundation Bearing Capacity Calculator

This script calculates the ultimate bearing capacity of shallow foundations
using Terzaghi's and Meyerhof's methods for different foundation shapes.

Author: Foundation Analysis Tool
"""

import math
import yaml
import json
import argparse
import sys
from pathlib import Path
from typing import Dict, Tuple, Optional, List, Union
from enum import Enum


class FoundationConfig:
    """Configuration class for foundation bearing capacity analysis."""
    
    def __init__(self, config_input: Union[str, Dict] = "config.yaml"):
        """
        Initialize configuration from YAML file, JSON file, or dictionary.
        
        Args:
            config_input: Path to the YAML/JSON configuration file or a dictionary
        """
        if isinstance(config_input, dict):
            # Direct dictionary input
            self.config_path = None
            self.config = config_input
        else:
            # File path input
            self.config_path = Path(config_input)
            self.config = self._load_config()
        self._validate_config()
    
    def _load_config(self) -> Dict:
        """Load configuration from YAML or JSON file."""
        if not self.config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {self.config_path}")
        
        file_extension = self.config_path.suffix.lower()
        
        try:
            with open(self.config_path, 'r') as file:
                if file_extension in ['.yaml', '.yml']:
                    config = yaml.safe_load(file)
                elif file_extension == '.json':
                    config = json.load(file)
                else:
                    # Try YAML first, then JSON as fallback
                    file.seek(0)
                    content = file.read()
                    try:
                        config = yaml.safe_load(content)
                    except yaml.YAMLError:
                        try:
                            config = json.loads(content)
                        except json.JSONDecodeError:
                            raise ValueError(f"Unable to parse configuration file as YAML or JSON: {self.config_path}")
                return config
        except (yaml.YAMLError, json.JSONDecodeError) as e:
            raise ValueError(f"Error parsing configuration file: {e}")
    
    def _validate_config(self) -> None:
        """Validate the loaded configuration."""
        required_sections = ['soil', 'foundation', 'loading', 'analysis', 'output']
        for section in required_sections:
            if section not in self.config:
                raise ValueError(f"Missing required section in config: {section}")
        
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
    def gwt(self) -> Optional[float]:
        return self.config['soil'].get('gwt')
    
    @property
    def foundation_width(self) -> float:
        return self.config['foundation']['width']
    
    @property
    def foundation_length(self) -> Optional[float]:
        return self.config['foundation'].get('length')
    
    @property
    def foundation_depth(self) -> float:
        return self.config['foundation']['depth']
    
    @property
    def load_angle(self) -> float:
        return self.config['loading'].get('load_angle', 0.0)
    
    @property
    def eccentricity(self) -> float:
        return self.config['loading'].get('eccentricity', 0.0)
    
    
    @property
    def safety_factor(self) -> float:
        return self.config['analysis']['safety_factor']
    
    @property
    def show_factors(self) -> bool:
        return self.config['output']['show_factors']
    
    @property
    def show_components(self) -> bool:
        return self.config['output']['show_components']
    
    @property
    def decimal_places(self) -> int:
        return self.config['output']['decimal_places']


class BearingCapacityCalculator:
    def __init__(self):
        self.bearing_capacity_factors = {}
    
    def _validate_positive_number(self, value: float, name: str) -> None:
        """Validate that a value is a positive number."""
        if not isinstance(value, (int, float)) or value < 0:
            raise ValueError(f"{name} must be a positive number, got {value}")
    
    def _validate_angle(self, angle: float, name: str) -> None:
        """Validate that an angle is within reasonable bounds."""
        if not isinstance(angle, (int, float)) or angle < 0 or angle > 50:
            raise ValueError(f"{name} must be between 0 and 50 degrees, got {angle}")
    
    def _validate_dimensions(self, width: float, length: float = None) -> None:
        """Validate foundation dimensions."""
        self._validate_positive_number(width, "Width")
        if length is not None:
            self._validate_positive_number(length, "Length")
            if width > length:
                raise ValueError("Width should not be greater than length for rectangular foundations")
    
    def get_bearing_capacity_factors(self, phi: float) -> Dict[str, float]:
        """
        Calculate bearing capacity factors Nc, Nq, and Nγ based on friction angle.
        
        Args:
            phi: Internal friction angle in degrees
            
        Returns:
            Dictionary containing Nc, Nq, and Nγ factors
        """
        self._validate_angle(phi, "Internal friction angle")
        phi_rad = math.radians(phi)
        
        # Calculate Nq
        Nq = math.exp(math.pi * math.tan(phi_rad)) * (math.tan(math.pi/4 + phi_rad/2))**2
        
        # Calculate Nc
        Nc = (Nq - 1) / math.tan(phi_rad) if phi > 0 else 5.7
        
        # Calculate Nγ (Meyerhof's formula)
        Nγ = (Nq - 1) * math.tan(1.4*phi_rad)
        
        return {
            'Nc': Nc,
            'Nq': Nq,
            'Nγ': Nγ
        }
    
    def get_shape_factors(self,
                         length: float, width: float ,
                         phi: float) -> Dict[str, float]:
        """
        Calculate shape factors for different foundation geometries.
        
        Args:
            length: Foundation length (for rectangular)
            width: Foundation width (for rectangular)
            phi: Internal friction angle in degrees
            
        Returns:
            Dictionary containing shape factors sc, sq, sγ
        """
        B_L = width / length
        Kp = math.tan(math.radians(45 + phi/2)) ** 2 # Passive earth pressure coefficient
        sc = 1.0 + 0.2 * Kp * B_L
        sq = 1.0 + 0.1 * Kp * B_L if phi > 10 else 1.0
        sγ = sq

        return {'sc': sc, 'sq': sq, 'sγ': sγ}
    
    def get_depth_factors(self, depth: float, width: float, phi: float) -> Dict[str, float]:
        """
        Calculate depth factors for foundations.
        
        Args:
            depth: Foundation depth
            width: Foundation width
            phi: Internal friction angle in degrees
            
        Returns:
            Dictionary containing depth factors dc, dq, dγ
        """
        D_B = depth / width
        Kp = math.tan(math.radians(45 + phi/2)) ** 2
        dc = 1 + 0.2 * math.sqrt(Kp) * D_B
        dq = 1 + 0.1 * math.sqrt(Kp) * D_B if phi > 10 else 1.0
        dγ = dq
        
        return {'dc': dc, 'dq': dq, 'dγ': dγ}

    def get_inclined_load_factors(self, load_angle: float, phi: float) -> Dict[str, float]:
        """
        Calculate factors for inclined loads.
        
        Args:
            load_angle: Load inclination angle in degrees
            phi: Internal friction angle in degrees
            
        Returns:
            Dictionary with load inclination factors ic, iq, iγ
        """
        self._validate_angle(load_angle, "Load angle")
        self._validate_angle(phi, "Internal friction angle")

        ic = (1 - load_angle / 90) ** 2
        iq = (1 - load_angle / 90) ** 2
        iγ = (1 - load_angle / phi) ** 2 if (phi > 10) else 1.0
        

        return {'ic': ic, 'iq': iq, 'iγ': iγ}
    
    def get_eccentricity_factors(self, e: float, width: float) -> Dict[str, float]:
        """
        Calculate factors for eccentric loads.

        Args:
            e (float): Eccentricity of load (m)
            width (float): Foundation width (m)

        Returns:
            Dict[str, float]: Dictionary with eccentricity factors
        """
        assert e >= 0, "Eccentricity must be non-negative"
        assert e < width/2, "Eccentricity must be less than half the width"
        ec = 1 - 2 * (e / width)
        eq = 1 - 2 * (e / width)
        eγ = 1 - 3.5 * (e / width) + 3 * (e / width) ** 2

        return {'ec': ec, 'eq': eq, 'eγ': eγ}
    
    def calculate_overburden_effective_stress(self, 
                                            unit_weights: List[float],
                                            layer_thicknesses: List[float],
                                            gwt_depth: Optional[float],
                                            foundation_depth: float) -> float:
        """
        Calculate overburden effective stress at foundation level for layered soil.
        
        Args:
            unit_weights: List of unit weights for each soil layer (kN/m³)
            layer_thicknesses: List of thicknesses for each soil layer (m)
            gwt_depth: Groundwater table depth from surface (m), None if no groundwater
            foundation_depth: Depth of foundation from surface (m)
            
        Returns:
            Overburden effective stress at foundation level (kPa)
        """
        if len(unit_weights) != len(layer_thicknesses):
            raise ValueError("Number of unit weights must match number of layer thicknesses")
        
        if foundation_depth < 0:
            raise ValueError("Foundation depth must be non-negative")
            
        total_stress = 0.0
        pore_pressure = 0.0
        current_depth = 0.0
        remaining_depth = foundation_depth
        
        # Water unit weight (kN/m³)
        gamma_w = 9.81
        
        for gamma, thickness in zip(unit_weights, layer_thicknesses):
            if remaining_depth <= 0:
                break
                
            # Determine the effective thickness of this layer
            layer_contribution = min(thickness, remaining_depth)
            
            # Calculate total stress contribution from this layer
            total_stress += gamma * layer_contribution
            
            # Calculate pore pressure if groundwater table is present
            if gwt_depth is not None:
                layer_top = current_depth
                layer_bottom = current_depth + layer_contribution
                
                if layer_bottom > gwt_depth:
                    # This layer is partially or fully below GWT
                    submerged_thickness = layer_bottom - max(layer_top, gwt_depth)
                    pore_pressure += gamma_w * submerged_thickness
            
            current_depth += layer_contribution
            remaining_depth -= layer_contribution
        
        # If foundation depth extends beyond all defined layers, use the last layer's unit weight
        if remaining_depth > 0:
            if not unit_weights:
                raise ValueError("No soil layers defined")
            
            last_gamma = unit_weights[-1]
            total_stress += last_gamma * remaining_depth
            
            # Continue pore pressure calculation if needed
            if gwt_depth is not None and current_depth >= gwt_depth:
                pore_pressure += gamma_w * remaining_depth
        
        # Effective stress = Total stress - Pore pressure
        effective_stress = total_stress - pore_pressure
        
        return max(0.0, effective_stress)  # Ensure non-negative result

    def calculate_effective_shear_depth(self,
                                      width: float,
                                      cohesions: List[float],
                                      friction_angles: List[float],
                                      layer_thicknesses: List[float],
                                      foundation_depth: float,
                                      tolerance: float = 0.0001,
                                      max_iterations: int = 100) -> Dict[str, float]:
        """
        Calculate effective shear depth (He) with iterative approach to find average soil parameters.
        
        He = 0.5 * B * tan(45 + φ/2)
        
        Args:
            width: Foundation width (m)
            cohesions: List of cohesion values for each layer (kPa)
            friction_angles: List of friction angles for each layer (degrees)
            layer_thicknesses: List of thicknesses for each layer (m)
            foundation_depth: Foundation depth from surface (m)
            tolerance: Convergence tolerance for He (default 0.0001)
            max_iterations: Maximum number of iterations (default 100)
            
        Returns:
            Dictionary containing He, cav, φav, and iteration count
        """
        if len(cohesions) != len(friction_angles) or len(cohesions) != len(layer_thicknesses):
            raise ValueError("Cohesions, friction angles, and thicknesses lists must have same length")
        
        if width <= 0 or foundation_depth < 0:
            raise ValueError("Width must be positive and foundation depth non-negative")
        
        # Find the layer just under the foundation
        current_depth = 0.0
        foundation_layer_index = 0
        
        for i, thickness in enumerate(layer_thicknesses):
            if current_depth + thickness > foundation_depth:
                foundation_layer_index = i
                break
            current_depth += thickness
        else:
            # Foundation is below all defined layers, use last layer
            foundation_layer_index = len(friction_angles) - 1
        
        # Initial He calculation using friction angle of layer under foundation
        phi_initial = friction_angles[foundation_layer_index]
        he_previous = 0.5 * width * math.tan(math.radians(45 + phi_initial / 2))
        
        iteration = 0
        converged = False
        
        while iteration < max_iterations and not converged:
            # Calculate average properties within current He
            cav, phi_av = self._calculate_average_properties(
                cohesions, friction_angles, layer_thicknesses,
                foundation_depth, he_previous
            )
            
            # Calculate new He using average friction angle
            he_new = 0.5 * width * math.tan(math.radians(45 + phi_av / 2))
            
            # Check convergence
            if abs(he_new - he_previous) < tolerance:
                converged = True
            
            he_previous = he_new
            iteration += 1
        
        if not converged:
            print(f"Warning: He calculation did not converge after {max_iterations} iterations")
        
        return {
            'He': he_previous,
            'cav': cav,
            'phi_av': phi_av,
            'iterations': iteration,
            'converged': converged
        }
    
    def _calculate_average_properties(self,
                                    cohesions: List[float],
                                    friction_angles: List[float],
                                    layer_thicknesses: List[float],
                                    foundation_depth: float,
                                    he: float) -> Tuple[float, float]:
        """
        Calculate weighted average cohesion and friction angle within effective shear depth.
        
        Args:
            cohesions: List of cohesion values (kPa)
            friction_angles: List of friction angles (degrees)
            layer_thicknesses: List of layer thicknesses (m)
            foundation_depth: Foundation depth (m)
            he: Effective shear depth (m)
            
        Returns:
            Tuple of (average cohesion, average friction angle)
        """
        total_weighted_c = 0.0
        total_weighted_phi = 0.0
        total_thickness = 0.0
        
        remaining_he = he
        
        # Start from foundation level and go down
        layer_start_depth = 0.0
        for i, thickness in enumerate(layer_thicknesses):
            layer_end_depth = layer_start_depth + thickness
            
            if layer_end_depth <= foundation_depth:
                # Layer is above foundation, skip
                layer_start_depth = layer_end_depth
                continue
            
            if remaining_he <= 0:
                # No more effective shear depth to consider
                break
            
            # Determine the thickness of this layer within the effective shear zone
            layer_contribution_start = max(layer_start_depth, foundation_depth)
            layer_contribution_end = min(layer_end_depth, foundation_depth + remaining_he)
            layer_contribution_thickness = max(0, layer_contribution_end - layer_contribution_start)
            
            if layer_contribution_thickness > 0:
                total_weighted_c += cohesions[i] * layer_contribution_thickness
                total_weighted_phi += math.tan(math.radians(friction_angles[i])) * layer_contribution_thickness
                total_thickness += layer_contribution_thickness
                remaining_he -= layer_contribution_thickness
            
            layer_start_depth = layer_end_depth
        
        # If He extends beyond all defined layers, use the last layer's properties
        if remaining_he > 0 and len(cohesions) > 0:
            total_weighted_c += cohesions[-1] * remaining_he
            total_weighted_phi += math.tan(math.radians(friction_angles[-1])) * remaining_he
            total_thickness += remaining_he
        
        if total_thickness == 0:
            raise ValueError("No soil layers found within effective shear depth")
        
        cav = total_weighted_c / total_thickness
        phi_av = math.degrees(math.atan(total_weighted_phi / total_thickness))

        return cav, phi_av

    def _calculate_average_unit_weight(self,
                                     unit_weights: List[float],
                                     layer_thicknesses: List[float],
                                     gwt_depth: Optional[float],
                                     foundation_depth: float,
                                     he: float) -> float:
        """
        Calculate area-weighted average effective unit weight for downward-tapering wedge.
        Uses the formula: γe = (1/H²) × Σ[γi × ((zi+1-zi)×H - 0.5×(zi+1² - zi²))]
        Uses effective unit weight (γ' = γ - γw) below groundwater table.
        
        Args:
            unit_weights: List of unit weights for each layer (kN/m³)
            layer_thicknesses: List of layer thicknesses (m)
            gwt_depth: Groundwater table depth from surface (m), None if no groundwater
            foundation_depth: Foundation depth from surface (m)
            he: Effective shear depth (m)
            
        Returns:
            Area-weighted average effective unit weight within effective shear depth wedge
        """
        if he <= 0:
            raise ValueError("Effective shear depth must be positive")
        
        # Water unit weight (kN/m³)
        gamma_w = 9.81
        
        # Build layer boundaries within the effective shear depth
        layer_boundaries = []  # Depths from foundation level (z values)
        current_depth = 0.0  # Depth from ground surface
        shear_zone_top = foundation_depth
        shear_zone_bottom = foundation_depth + he
        # Add foundation level as starting point
        if foundation_depth >= current_depth:
            layer_boundaries.append(0.0)  # z = 0 at foundation level
        
        # Process each soil layer
        for i, thickness in enumerate(layer_thicknesses):
            layer_top = current_depth
            layer_bottom = current_depth + thickness
            
            # Check if this layer intersects with the effective shear depth zone
            if layer_bottom > shear_zone_top and layer_top < shear_zone_bottom:
                # Layer intersects with shear zone
                intersection_top = max(layer_top, shear_zone_top)
                intersection_bottom = min(layer_bottom, shear_zone_bottom)
                
                # Convert to depths from foundation level
                z_top = intersection_top - foundation_depth
                z_bottom = intersection_bottom - foundation_depth
                
                # Add boundaries if not already added
                if z_top not in layer_boundaries:
                    layer_boundaries.append(z_top)
                if z_bottom not in layer_boundaries:
                    layer_boundaries.append(z_bottom)
            
            current_depth += thickness
        
        # Ensure we have the full effective shear depth
        if he not in layer_boundaries:
            layer_boundaries.append(he)
        
        # Sort boundaries
        layer_boundaries.sort()
        # Calculate weighted sum using the wedge formula
        weighted_sum = 0.0
        
        for i in range(len(layer_boundaries) - 1):
            zi = layer_boundaries[i]
            zi_plus_1 = layer_boundaries[i + 1]
            
            if zi >= he:  # Beyond effective shear depth
                break
            
            # Find which soil layer this segment belongs to
            segment_mid_abs = foundation_depth + (zi + zi_plus_1) / 2
            current_layer_depth = 0.0
            gamma_layer = unit_weights[0]  # Default to first layer
            
            for j, thickness in enumerate(layer_thicknesses):
                if current_layer_depth + thickness > segment_mid_abs:
                    gamma_layer = unit_weights[j]
                    break
                current_layer_depth += thickness
            else:
                # If beyond all defined layers, use last layer
                if len(unit_weights) > 0:
                    gamma_layer = unit_weights[-1]
            
            # Apply the wedge formula for this layer segment
            # γi × [(zi+1-zi)×H - 0.5×(zi+1² - zi²)]
            layer_contribution = gamma_layer * (
                (zi_plus_1 - zi) * he - 0.5 * (zi_plus_1**2 - zi**2)
            )
            weighted_sum += layer_contribution

        # Calculate dw (GWT - foundation depth >= 0)
        dw = (gwt_depth - foundation_depth) if gwt_depth is not None else he
        # Ensure non-negative effective dw
        dw = max(0.0, dw)
        # Calculate the water contribution
        water_contribution = gamma_w * (
                (he - dw) * he - 0.5 * (he**2 - dw**2)
            ) if he-dw > 0 else 0.0
        # Apply the final formula: γe = (1/H²) × weighted_sum
        weighted_sum -= water_contribution
        gamma_average = 2 * weighted_sum / (he**2)
        
        return max(0.0, gamma_average)  # Ensure non-negative result

    def print_soil_profile_diagram(self,
                                  unit_weights: List[float],
                                  cohesions: List[float],
                                  friction_angles: List[float],
                                  layer_thicknesses: List[float],
                                  gwt_depth: Optional[float],
                                  foundation_width: float,
                                  foundation_depth: float,
                                  he: float) -> None:
        """
        Print a schematic ASCII diagram of the soil profile, foundation, and effective shear depth.
        Shows only important depths (not to scale) for clarity.
        
        Args:
            unit_weights: List of unit weights for each layer (kN/m³)
            cohesions: List of cohesion values for each layer (kPa)
            friction_angles: List of friction angles for each layer (degrees)
            layer_thicknesses: List of thicknesses for each layer (m)
            gwt_depth: Groundwater table depth from surface (m)
            foundation_width: Foundation width (m)
            foundation_depth: Foundation depth (m)
            he: Effective shear depth (m)
        """
        print("\n" + "="*60)
        print("SOIL PROFILE DIAGRAM (Schematic - Not to Scale)")
        print("="*60)
        
        # Identify important depths
        important_depths = [0.0]  # Ground surface
        
        # Add foundation depth
        important_depths.append(foundation_depth)
        
        # Add GWT if present
        if gwt_depth is not None:
            important_depths.append(gwt_depth)
        
        # Add layer boundaries within relevant zone
        current_depth = 0.0
        analysis_zone = foundation_depth + he + 1.0  # Show a bit beyond He
        
        for thickness in layer_thicknesses:
            current_depth += thickness
            if current_depth <= analysis_zone:
                important_depths.append(current_depth)
        
        # Add He bottom
        he_bottom = foundation_depth + he
        important_depths.append(he_bottom)
        
        # Remove duplicates and sort
        important_depths = sorted(list(set(important_depths)))
        
        diagram_width = 50
        
        print(f"{'Depth':>8} | {'Profile':^{diagram_width}} |")
        print(f"{'-'*8}-+-{'-'*diagram_width}-+")
        
        for i, depth in enumerate(important_depths):
            # Determine what's at this depth
            layer_info = self._get_layer_at_depth(depth, layer_thicknesses, unit_weights, cohesions, friction_angles)
            layer_symbol = layer_info['symbol']
            
            # Create the line representation
            line = [' '] * diagram_width
            
            # Ground surface
            if depth == 0.0:
                for j in range(diagram_width):
                    line[j] = '='
            
            # Layer boundaries
            elif i > 0 and depth in [sum(layer_thicknesses[:k+1]) for k in range(len(layer_thicknesses))]:
                for j in range(diagram_width):
                    line[j] = '-'
                # Add layer symbol
                line[2] = layer_symbol
            
            # Foundation level
            elif abs(depth - foundation_depth) < 0.01:  # Foundation depth
                found_center = diagram_width // 2
                found_half_width = min(8, diagram_width // 4)
                
                # Check if GWT coincides with foundation depth
                if gwt_depth is not None and abs(depth - gwt_depth) < 0.01:
                    # Show both foundation and GWT
                    for j in range(diagram_width):
                        line[j] = '~'
                    # Foundation base over GWT
                    for j in range(found_center - found_half_width, found_center + found_half_width + 1):
                        if 0 <= j < diagram_width:
                            line[j] = '#'
                else:
                    # Foundation base only
                    for j in range(found_center - found_half_width, found_center + found_half_width + 1):
                        if 0 <= j < diagram_width:
                            line[j] = '#'
                
                # Add layer symbol
                line[2] = layer_symbol
            
            # Groundwater table
            elif gwt_depth is not None and abs(depth - gwt_depth) < 0.01:
                for j in range(diagram_width):
                    line[j] = '~'
                # Add layer symbol
                line[2] = layer_symbol
            
            # He bottom (show wedge bottom)
            elif abs(depth - he_bottom) < 0.01:
                found_center = diagram_width // 2
                wedge_half_width = min(12, diagram_width // 3)
                
                # Wedge outline
                for j in [found_center - wedge_half_width, found_center + wedge_half_width]:
                    if 0 <= j < diagram_width:
                        line[j] = '.'
                
                # Add layer symbol
                line[2] = layer_symbol
            
            # Regular layer line
            else:
                line[2] = layer_symbol
            
            # Add foundation sides if within foundation zone
            if foundation_depth <= depth <= he_bottom:
                found_center = diagram_width // 2
                found_half_width = min(8, diagram_width // 4)
                
                # Add wedge outline if within He zone
                if depth > foundation_depth:
                    # Calculate wedge width based on depth
                    depth_from_found = depth - foundation_depth
                    wedge_ratio = min(1.0, depth_from_found / he)
                    wedge_half_width = int(wedge_ratio * min(12, diagram_width // 3))
                    
                    for j in [found_center - wedge_half_width, found_center + wedge_half_width]:
                        if 0 <= j < diagram_width and line[j] == ' ':
                            line[j] = '.'
            
            line_str = ''.join(line)
            print(f"{depth:7.2f}m | {line_str} |")
        
        # Print legend
        print(f"{'-'*8}-+-{'-'*diagram_width}-+")
        print("\nLEGEND:")
        print("-"*40)
        print("=== : Ground Surface")
        print("--- : Layer Boundary")
        print("~~~ : Groundwater Table")
        print("### : Foundation")
        print("... : Effective Shear Depth Wedge")
        print("")
        
        # Print layer properties
        print("LAYER PROPERTIES:")
        print("-"*60)
        current_depth = 0.0
        for i, (thickness, gamma, c, phi) in enumerate(zip(layer_thicknesses, unit_weights, cohesions, friction_angles)):
            layer_symbol = chr(ord('A') + i) if i < 26 else str(i)
            depth_range = f"{current_depth:.1f}-{current_depth + thickness:.1f}m"
            print(f"Layer {layer_symbol} ({depth_range}): γ={gamma:.1f} kN/m³, c={c:.1f} kPa, φ={phi:.1f}°")
            current_depth += thickness
        
        print("")
        print("FOUNDATION & ANALYSIS:")
        print("-"*60)
        print(f"Foundation Depth: {foundation_depth:.1f} m")
        print(f"Foundation Width: {foundation_width:.1f} m")
        print(f"Effective Shear Depth (He): {he:.2f} m")
        if gwt_depth is not None:
            print(f"Groundwater Table: {gwt_depth:.1f} m")
        else:
            print("Groundwater Table: Not present")
        print("="*60)
    
    def _get_layer_at_depth(self, depth: float, layer_thicknesses: List[float], 
                           unit_weights: List[float], cohesions: List[float], 
                           friction_angles: List[float]) -> dict:
        """Get layer information at a specific depth."""
        current_depth = 0.0
        for i, thickness in enumerate(layer_thicknesses):
            if current_depth + thickness > depth:
                return {
                    'index': i,
                    'symbol': chr(ord('A') + i) if i < 26 else str(i),
                    'unit_weight': unit_weights[i],
                    'cohesion': cohesions[i],
                    'friction_angle': friction_angles[i]
                }
            current_depth += thickness
        
        # If beyond all layers, return last layer
        last_idx = len(layer_thicknesses) - 1
        return {
            'index': last_idx,
            'symbol': chr(ord('A') + last_idx) if last_idx < 26 else str(last_idx),
            'unit_weight': unit_weights[last_idx] if unit_weights else 0,
            'cohesion': cohesions[last_idx] if cohesions else 0,
            'friction_angle': friction_angles[last_idx] if friction_angles else 0
        }

    def calculate_bearing_capacity(self,
                                         cohesion: List[float],
                                         unit_weight: List[float],
                                         friction_angle: List[float],
                                         thickness: List[float],
                                         gwt_depth: Optional[float],
                                         width: float,
                                         depth: float,
                                         load_angle: float,
                                         eccentricity: float,
                                         length: float = None) -> Dict[str, float]:
        """
        Calculate ultimate bearing capacity using Meyerhof's method.
        
        Args:
            cohesion: Soil cohesion (kPa)
            unit_weight: Soil unit weight (kN/m³)
            friction_angle: Internal friction angle (degrees)
            width: Foundation width (m)
            depth: Foundation depth (m)
            length: Foundation length (m, for rectangular only)
            
        Returns:
            Dictionary with bearing capacity components and total
        """
        # Input validation
        if isinstance(cohesion, list):
            for c in cohesion:
                self._validate_positive_number(c, "Cohesion")
        else:
            self._validate_positive_number(cohesion, "Cohesion")
            
        if isinstance(unit_weight, list):
            for uw in unit_weight:
                self._validate_positive_number(uw, "Unit weight")
        else:
            self._validate_positive_number(unit_weight, "Unit weight")
            
        if isinstance(friction_angle, list):
            for fa in friction_angle:
                self._validate_angle(fa, "Friction angle")
        else:
            self._validate_angle(friction_angle, "Friction angle")
            
        self._validate_positive_number(depth, "Depth")
        self._validate_dimensions(width, length)
        
        He_result = self.calculate_effective_shear_depth(
            width=width,
            cohesions=cohesion,
            friction_angles=friction_angle,
            layer_thicknesses=thickness,
            foundation_depth=depth
        )
        
        # Use average cohesion and friction angle from effective shear depth
        cohesion = He_result['cav']
        friction_angle = He_result['phi_av']
        print(f"Effective shear depth (He): {He_result['He']:.2f} m")
        print(f"Average cohesion (cav): {cohesion:.2f} kPa")
        print(f"Average friction angle (φav): {friction_angle:.2f} degrees")
        print(f"Iterations: {He_result['iterations']}, Converged: {He_result['converged']}")

        # Get bearing capacity factors
        factors = self.get_bearing_capacity_factors(friction_angle)
        Nc, Nq, Nγ = factors['Nc'], factors['Nq'], factors['Nγ']
        
        # Get shape factors
        shape_factors = self.get_shape_factors(length, width, friction_angle)
        sc, sq, sγ = shape_factors['sc'], shape_factors['sq'], shape_factors['sγ']
        
        # Get depth factors
        depth_factors = self.get_depth_factors(depth, width, friction_angle)
        dc, dq, dγ = depth_factors['dc'], depth_factors['dq'], depth_factors['dγ']
        
        # Get inclined load factors
        inclined_factors = self.get_inclined_load_factors(load_angle, friction_angle)
        ic, iq, iγ = inclined_factors['ic'], inclined_factors['iq'], inclined_factors['iγ']
        if load_angle != 0.0:
            sc = 1.0
            sq = 1.0
            sγ = 1.0
        # Get eccentricity factors
        eccentricity_factors = self.get_eccentricity_factors(eccentricity, width)
        ec, eq, eγ = eccentricity_factors['ec'], eccentricity_factors['eq'], eccentricity_factors['eγ']
        
        # Get overburden pressure
        overburden_pressure = self.calculate_overburden_effective_stress(
            unit_weights=unit_weight,
            layer_thicknesses=thickness,
            gwt_depth=gwt_depth,
            foundation_depth=depth
        )
        
        # Get average unit weight within effective shear depth
        avg_unit_weight = self._calculate_average_unit_weight(
            unit_weights=unit_weight,
            layer_thicknesses=thickness,
            gwt_depth=gwt_depth,
            foundation_depth=depth,
            he=He_result['He']
        )

        # Meyerhof's bearing capacity equation
        qc = cohesion * Nc * sc * dc * ic * ec
        qq = overburden_pressure * Nq * sq * dq * iq * eq
        qγ = 0.5 * avg_unit_weight * width * Nγ * sγ * dγ * iγ * eγ
        print(f"qc = {cohesion} * {Nc:.2f} * {sc:.2f} * {dc:.2f} * {ic:.2f} * {ec:.2f} = {qc:.2f} kPa")
        print(f"qq = {overburden_pressure} * {Nq:.2f} * {sq:.2f} * {dq:.2f} * {iq:.2f} * {eq:.2f} = {qq:.2f} kPa")
        print(f"qγ = 0.5 * {avg_unit_weight:.2f} * {width} * {Nγ:.2f} * {sγ:.2f} * {dγ:.2f} * {iγ:.2f} * {eγ:.2f} = {qγ:.2f} kPa")
        qu = qc + qq + qγ
        
        return {
            'cohesion_term': qc,
            'surcharge_term': qq,
            'weight_term': qγ,
            'ultimate_bearing_capacity': qu,
            'bearing_capacity_factors': factors,
            'shape_factors': shape_factors,
            'depth_factors': depth_factors,
            'effective_shear_depth': He_result['He']
        }
    
    def calculate_allowable_bearing_capacity(self, ultimate_capacity: float, 
                                           safety_factor: float = 3.0) -> float:
        """
        Calculate allowable bearing capacity using factor of safety.
        
        Args:
            ultimate_capacity: Ultimate bearing capacity
            safety_factor: Factor of safety (default 3.0)
            
        Returns:
            Allowable bearing capacity
        """
        self._validate_positive_number(ultimate_capacity, "Ultimate capacity")
        if safety_factor < 1.0:
            raise ValueError(f"Safety factor must be >= 1.0, got {safety_factor}")
        
        return ultimate_capacity / safety_factor


def run_analysis_from_config(config_path: str) -> None:
    """
    Run bearing capacity analysis using configuration file.
    
    Args:
        config_path: Path to the YAML configuration file
    """
    try:
        # Load configuration
        config = FoundationConfig(config_path)
        calculator = BearingCapacityCalculator()
        
        print("Foundation Bearing Capacity Calculator")
        print("=" * 50)
        print(f"Configuration file: {config.config_path}")
        print("=" * 50)
        
        
        results = {}

        # Run Meyerhof method
        meyerhof_result = calculator.calculate_bearing_capacity(
            cohesion=config.soil_cohesion,
            unit_weight=config.soil_unit_weight,
            friction_angle=config.soil_friction_angle,
            thickness=config.soil_thickness,
            gwt_depth=config.gwt,
            width=config.foundation_width,
            depth=config.foundation_depth,
            load_angle=config.load_angle,
            eccentricity=config.eccentricity,
            length=config.foundation_length
        )
        results['meyerhof'] = meyerhof_result
        
        # Print soil profile diagram
        calculator.print_soil_profile_diagram(
            unit_weights=config.soil_unit_weight,
            cohesions=config.soil_cohesion,
            friction_angles=config.soil_friction_angle,
            layer_thicknesses=config.soil_thickness,
            gwt_depth=config.gwt,
            foundation_width=config.foundation_width,
            foundation_depth=config.foundation_depth,
            he=meyerhof_result['effective_shear_depth']
        )
        
        if config.show_components:
            print(f"   Cohesion term: {meyerhof_result['cohesion_term']:.{config.decimal_places}f} kPa")
            print(f"   Surcharge term: {meyerhof_result['surcharge_term']:.{config.decimal_places}f} kPa")
            print(f"   Weight term: {meyerhof_result['weight_term']:.{config.decimal_places}f} kPa")
        
        print(f"   Ultimate bearing capacity: {meyerhof_result['ultimate_bearing_capacity']:.{config.decimal_places}f} kPa")
        
        allowable_meyerhof = calculator.calculate_allowable_bearing_capacity(
            meyerhof_result['ultimate_bearing_capacity'], config.safety_factor
        )
        print(f"   Allowable bearing capacity (FS={config.safety_factor}): {allowable_meyerhof:.{config.decimal_places}f} kPa")
        
        # Display bearing capacity factors (using average friction angle from He calculation)
        if config.show_factors:
            # Get the average friction angle from the result
            avg_phi = meyerhof_result['bearing_capacity_factors']
            print(f"\nBearing Capacity Factors:")
            print(f"   Nc = {avg_phi['Nc']:.{config.decimal_places}f}")
            print(f"   Nq = {avg_phi['Nq']:.{config.decimal_places}f}")
            print(f"   Nγ = {avg_phi['Nγ']:.{config.decimal_places}f}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Foundation Bearing Capacity Calculator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python bearing_capacity.py                     # Use default config.yaml
  python bearing_capacity.py -c my_config.yaml  # Use custom YAML config file
  python bearing_capacity.py -c my_config.json  # Use custom JSON config file
        """
    )
    
    parser.add_argument(
        '-c', '--config',
        default='config.yaml',
        help='Path to YAML or JSON configuration file (default: config.yaml)'
    )
    
    args = parser.parse_args()
    
    # Run analysis with configuration
    run_analysis_from_config(args.config)


if __name__ == "__main__":
    main()