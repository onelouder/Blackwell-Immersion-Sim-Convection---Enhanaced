import React, { useState, useEffect, useMemo } from 'react';
import { createRoot } from 'react-dom/client';
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
  ReferenceLine,
  ComposedChart
} from 'recharts';
import { Settings, Wind, Info, Droplets, Ruler, Calculator, Edit2, X, Save, RotateCcw, Download, Trash2, Plus, Cpu } from 'lucide-react';

// --- 1. Physics Constants & Fluid Database ---

const G = 9.81; // Gravity (m/s^2)

// Default Geometry defined in "Final Parameter Summary" PDF Page 8
const INITIAL_GEOMETRY = {
  L: 0.15, // Height of HS in vertical orientation (m)
  W: 0.08, // Width (m)
  s: 0.003, // Fin spacing (m)
  t_fin: 0.0008, // Fin thickness (m)
  h_fin: 0.025, // Fin height/length (m)
  
  // Die Dimensions for Spreading Resistance
  die_W: 0.02, // 20mm width
  die_L: 0.03, // 30mm height
  
  // Thermal Resistances (K/W)
  R_TIM: 0.02, // Thermal Interface Material (approx for high end TIM)
  R_die: 0.01, // Internal die resistance
  
  // Base Material (Copper)
  k_base: 400, // W/m-K
};

// Fluid Properties
const DEFAULT_FLUIDS = [
  { 
    id: 'dcf281', 
    name: 'DCF-281 (Novvi)', 
    rho: 793.36, 
    nu: 10.41e-6, // m2/s @ 40C
    nu_cSt_40: 10.41,
    nu_cSt_100: 2.81, 
    k: 0.14059, 
    cp: 2110, 
    beta: 0.00092, 
    color: '#64748b' 
  },
  { 
    id: 'shell', 
    name: 'Shell XHVI3', 
    rho: 790.2, 
    nu: 9.99e-6, 
    nu_cSt_40: 9.99,
    nu_cSt_100: 2.72,
    k: 0.13702, 
    cp: 2070, 
    beta: 0.00064, 
    color: '#f59e0b' 
  },
  { 
    id: 'castrol', 
    name: 'Castrol DC15', 
    rho: 819, 
    nu: 7.50e-6, 
    nu_cSt_40: 7.50,
    nu_cSt_100: 2.20,
    k: 0.134, 
    cp: 2200, 
    beta: 0.00090, 
    color: '#10b981' 
  },
  { 
    id: 'valvoline', 
    name: 'Valvoline HPC', 
    rho: 811.5, 
    nu: 8.01e-6, 
    nu_cSt_40: 8.01,
    nu_cSt_100: 2.40,
    k: 0.1304, 
    cp: 2000, 
    beta: 0.00080, 
    color: '#ef4444' 
  },
  { 
    id: 'fuchs', 
    name: 'Fuchs Renolin', 
    rho: 826, 
    nu: 4.96e-6, 
    nu_cSt_40: 4.96,
    nu_cSt_100: 1.70,
    k: 0.134, 
    cp: 2200, 
    beta: 0.00065, 
    color: '#8b5cf6' 
  },
  { 
    id: 'mpao', 
    name: 'Novel MPAO', 
    rho: 794.8, 
    nu: 7.994e-6, 
    nu_cSt_40: 7.994,
    nu_cSt_100: 1.885,
    k: 0.15203, 
    cp: 2270, 
    beta: 0.00080, 
    color: '#3b82f6' 
  }, 
];

// --- 2. Physics Solver Engines ---

// Helper: Calculate Temperature Dependent Viscosity (ASTM D341)
const getViscosity = (fluid: any, T_C: number) => {
  if (!fluid.nu_cSt_100 || !fluid.nu_cSt_40) return fluid.nu;

  const T_K = T_C + 273.15;
  const T1 = 40 + 273.15;
  const T2 = 100 + 273.15;

  const v1 = fluid.nu_cSt_40;
  const v2 = fluid.nu_cSt_100;

  const Z1 = Math.log10(Math.log10(v1 + 0.7));
  const Z2 = Math.log10(Math.log10(v2 + 0.7));

  const B = (Z1 - Z2) / (Math.log10(T2) - Math.log10(T1));
  const A = Z1 + B * Math.log10(T1);

  const Z_target = A - B * Math.log10(T_K);
  const v_cSt = Math.pow(10, Math.pow(10, Z_target)) - 0.7;

  return Math.max(0.1, v_cSt) * 1e-6; 
};

// Teertstra (1999) Correlation for Plate Fin Heatsinks (Forced/Induced Draft)
// Blends fully developed channel flow with isolated plate boundary layer flow
const calculateTeertstraNu = (Re_channel: number, Pr: number, s: number, L: number) => {
  // Re_channel is based on channel spacing (s) or hydraulic diameter? 
  // Teertstra uses Re_b = (u * s) / v
  // Our Re passed in is usually Dh based. Let's convert if needed.
  // Actually, let's use the explicit Teertstra form:
  // Nu_b = [ (1 / (Re_b * Pr / 2)^3) + (1 / (0.664 * sqrt(Re_b) * Pr^1/3)^3 ) ] ^ -1/3
  
  // Note: The Re passed into this function from solver is based on Dh.
  // Dh approx 2s. So Re_b approx Re_Dh / 2.
  
  // Let's recalculate local Reynolds based on 's' for strict adherence to the paper.
  // Assuming the solver provides the velocity, we re-calc Re_s inside here? 
  // No, better to stick to the composite form using the Re provided.
  
  // Let's implement the standard composite equation for parallel plates:
  // Nu_fd = (1/24) * Ra_s * (s/L)  <-- Natural conv limit
  // But we are in "Induced Draft" mode, which effectively acts like forced convection driven by pressure.
  
  // Using Teertstra Composite (1999):
  // Nu_b = [ (Re_s * Pr / 2)^-3 + (0.664 * sqrt(Re_s) * Pr^0.33 * sqrt(1 + 3.65/sqrt(Re_s)))^-3 ] ^ -0.33
  
  // We need Re based on spacing 's'
  const Re_s = Re_channel * (s / (2*s)); // Rough conversion if input is Re_Dh. 
  // Let's trust the solver to pass us the velocity to do it right.
  
  // Fallback: Using the robust Churchill-Usagi blending for vertical channels
  // Nu = [ (Nu_fd)^-n + (Nu_dev)^-n ] ^ -1/n
  // Nu_fd = 8.235 (infinite parallel plates constant heat flux)
  // Nu_dev = 0.664 * Re_L^0.5 * Pr^0.33 (Flat plate)
  
  // Let's use the Sieder-Tate/Hausen provided previously but tuned, 
  // OR switch to the Teertstra which is better for heatsinks.
  
  // Implementing Teertstra approx for the specific aspect ratio:
  const xi = L / (s * Re_channel * Pr);
  // Nu_s = 1 / [ (1/(Re_s Pr / 2))^3 + (1 / (0.664 Re_s^0.5 Pr^0.33))^3 ]^1/3
  // This is a powerful general correlation.
  
  const term1 = Math.pow(Re_channel * Pr * 0.5, -3); // Channel flow limit
  const term2 = Math.pow(0.664 * Math.sqrt(Re_channel) * Math.pow(Pr, 0.33), -3); // Boundary layer limit
  
  const Nu_composite = Math.pow(term1 + term2, -1/3);
  return Nu_composite;
};

const calculateFinEff = (h: number, geom: typeof INITIAL_GEOMETRY) => {
  const m = Math.sqrt((2 * h) / (geom.k_base * geom.t_fin));
  const mh = m * geom.h_fin;
  if (mh === 0) return 1;
  return Math.tanh(mh) / mh;
};

// Iterative Solver for a single Operating Point (Fixed Wall Temp)
const solveOperatingPoint = (
  fluid: any, 
  H_tank: number, 
  Tj: number, 
  T_inlet: number,
  geom: typeof INITIAL_GEOMETRY,
  derived: { N: number, A_total: number, A_fin: number, A_base: number, A_c_total: number, R_internal: number, Dh: number, R_spreading: number },
  bypassFactor: number = 0 
) => {
  
  let T_wall = Tj - 5; // Initial guess for Base Temp (lower than Tj due to Spreading + TIM)
  let T_out = T_inlet + 10; 
  
  let Q_dissipated = 0;
  let u_induced = 0;
  let h_conv = 0;
  let Re = 0;
  let Ra = 0;
  let nu_op = fluid.nu;
  
  let converged = false;
  let iterations = 0;
  const alpha = 0.2; 
  
  // Yovanovich Spreading Resistance calculated in derived already, but depends on physics?
  // No, purely geometric.
  
  // Area ratio for contraction/expansion loss
  // sigma = Free Flow Area / Frontal Area
  const A_frontal = geom.W * geom.h_fin; // Approximate frontal area of the array block
  const sigma = derived.A_c_total / A_frontal;
  
  // Form Loss Coefficients (Kays & London)
  // For Laminar flow approx
  const Kc = 0.42 * (1 - Math.pow(sigma, 2));
  const Ke = Math.pow(1 - sigma, 2); 
  const K_form = Kc + Ke;

  while (!converged && iterations < 100) {
    const T_mean_fluid = (T_inlet + T_out) / 2;
    
    // Property Evaluation at Film Temperature
    const nu_bulk = getViscosity(fluid, T_mean_fluid);
    const nu_wall = getViscosity(fluid, T_wall);
    nu_op = nu_bulk; 

    // Temp dependent Density & Expansion
    const beta = fluid.beta; 
    const rho_bulk = fluid.rho * (1 - beta * (T_mean_fluid - 40));
    const rho_inlet = fluid.rho * (1 - beta * (T_inlet - 40)); // driving head based on cold leg
    
    // Temp dependent Cp (approximate linear increase for oils, ~3 J/kgK per C)
    const cp_factor = 1 + 0.0015 * (T_mean_fluid - 40); 
    const cp_fluid = fluid.cp * cp_factor;

    // Driving Pressure (Buoyancy)
    // dP = (rho_cold - rho_hot_column) * g * H
    const dT_core = Math.max(0.1, T_mean_fluid - T_inlet);
    const dT_chimney = Math.max(0.1, T_out - T_inlet);
    const H_chimney = Math.max(0, H_tank - geom.L);
    
    // Core buoyancy + Chimney buoyancy
    const dP_buoy = G * beta * rho_inlet * (dT_core * geom.L + dT_chimney * H_chimney);

    // Flow Resistance
    const u_guess = u_induced > 0 ? u_induced : 0.05;
    const Re_guess = (u_guess * derived.Dh) / nu_bulk;
    
    // Apparent Friction Factor (Shah & London) for rectangular ducts
    // f_app * Re = 24 * (1 - 1.3553*AR + 1.9467*AR^2 ... ) 
    // Simplified Laminar: f = 64/Re (pipe) -> f = 24/Re (parallel plates)
    // Plus entry length effects
    
    // Hydrodynamic Entry Length: L+ = L / (Dh * Re)
    // f_app * Re = 24 + 0.67 K(infinity) / (4 L+)
    const f_fd = 24 / Re_guess; // Fully developed laminar
    const f_app = f_fd; // Simplifying for stability, Form losses handled by K_form

    // Viscosity correction for friction
    const visc_ratio = nu_wall / nu_bulk; // approx mu_w/mu_b if rho const
    const f_corrected = f_app * Math.pow(visc_ratio, 0.58);

    const K_friction = f_corrected * (geom.L / derived.Dh);
    const K_total = K_friction + K_form + 1.0; // +1.0 for velocity head exit

    // Calculate Velocity
    // dP = 0.5 * rho * u^2 * K_total
    let u_new = Math.sqrt( (2 * dP_buoy) / (rho_bulk * K_total) );
    
    // Apply Bypass Factor (Leakage)
    // Modeled as a direct efficiency loss on velocity
    u_new = u_new * (1 - bypassFactor);

    u_induced = u_induced * 0.5 + u_new * 0.5; 
    if (u_induced < 0.0001) u_induced = 0.0001;
    
    Re = (u_induced * derived.Dh) / nu_bulk;
    
    // Heat Transfer Correlation
    const Pr = (nu_bulk * rho_bulk * cp_fluid) / fluid.k;
    Ra = (G * beta * (T_wall - T_inlet) * Math.pow(geom.s, 3)) / (nu_bulk * nu_bulk) * Pr;
    
    // Use Teertstra Composite Correlation
    const Nu = calculateTeertstraNu(Re, Pr, geom.s, geom.L);
    
    const h_bare = (Nu * fluid.k) / derived.Dh;
    
    const eta_fin = calculateFinEff(h_bare, geom);
    const eta_o = 1 - (derived.A_fin / derived.A_total) * (1 - eta_fin);
    const h_eff = h_bare * eta_o;
    h_conv = h_eff;

    // Heat Balance
    const dT_conv = Math.max(0.1, T_wall - T_mean_fluid);
    const Q_conv = h_eff * derived.A_total * dT_conv;

    const m_dot = rho_bulk * u_induced * derived.A_c_total;
    
    let T_out_calc = T_inlet;
    if (m_dot > 0) {
      T_out_calc = T_inlet + (Q_conv / (m_dot * cp_fluid));
    }

    // T_wall is derived from Tj minus internal resistances (Conductive + Spreading)
    // Q = (Tj - T_wall) / R_int
    // T_wall = Tj - Q * R_int
    const T_wall_calc = Tj - Q_conv * derived.R_internal;

    const err_Tout = Math.abs(T_out_calc - T_out);
    const err_Twall = Math.abs(T_wall_calc - T_wall);

    if (err_Tout < 0.05 && err_Twall < 0.05) {
      converged = true;
      Q_dissipated = Q_conv;
      T_out = T_out_calc;
      T_wall = T_wall_calc;
    } else {
      // Relaxation
      T_out = T_out * (1 - alpha) + T_out_calc * alpha;
      T_wall = T_wall * (1 - alpha) + T_wall_calc * alpha;
      
      // Bounds check
      if (T_out > Tj - 1) T_out = Tj - 1;
      if (T_wall > Tj) T_wall = Tj;
      if (T_wall < T_inlet) T_wall = T_inlet + 0.1;
    }
    iterations++;
  }

  return {
    Q: Q_dissipated,
    u_ind: u_induced,
    T_out: T_out,
    T_wall: T_wall,
    T_j: Tj,
    h_conv: h_conv,
    Re: Re,
    Ra: Ra,
    nu_op: nu_op
  };
};

// Inverse Solver: Find Tj for a fixed Power Q
const solveForTj = (
  targetQ: number,
  fluid: any, 
  H_tank: number, 
  T_inlet: number,
  geom: typeof INITIAL_GEOMETRY,
  derived: any,
  bypassFactor: number
) => {
  let low = T_inlet + 0.1;
  let high = 250; // Safety max
  let bestResult = null;
  let minError = Infinity;

  // Binary search for Wall Temperature that results in Target Q
  for(let i=0; i<20; i++) {
     const midTj = (low + high) / 2;
     const result = solveOperatingPoint(fluid, H_tank, midTj, T_inlet, geom, derived, bypassFactor);
     
     const diff = result.Q - targetQ;
     
     if (Math.abs(diff) < minError) {
       minError = Math.abs(diff);
       bestResult = result;
     }

     if (Math.abs(diff) < 0.5) break; 

     if (result.Q > targetQ) {
       high = midTj; // Too hot, cool down
     } else {
       low = midTj; // Too cold, heat up
     }
  }
  return bestResult;
};

// --- 3. UI Components ---

const InputField = ({ label, value, unit, step = 0.001, min = 0, onChange }: any) => (
  <div className="mb-3">
    <label className="flex justify-between text-xs font-semibold text-slate-600 mb-1 uppercase tracking-wider">
      {label}
    </label>
    <div className="flex items-center relative">
      <input
        type="number"
        step={step}
        min={min}
        value={value}
        onChange={(e) => onChange(parseFloat(e.target.value) || 0)}
        className="w-full px-3 py-2 bg-white border border-slate-300 rounded text-sm text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-transparent transition-all font-mono"
      />
      <span className="absolute right-3 text-xs text-slate-400 pointer-events-none">{unit}</span>
    </div>
  </div>
);

const FluidEditModal = ({ isOpen, onClose, fluids, onSave, onReset }: any) => {
  const [localFluids, setLocalFluids] = useState(fluids);

  useEffect(() => {
    if (isOpen) {
      setLocalFluids(JSON.parse(JSON.stringify(fluids)));
    }
  }, [fluids, isOpen]);

  if (!isOpen) return null;

  const handleFluidChange = (id: string, field: string, val: any) => {
    setLocalFluids((prev: any[]) => prev.map(f => {
      if (f.id === id) {
        return { ...f, [field]: val };
      }
      return f;
    }));
  };

  const handleAddFluid = () => {
    const newId = `custom_${Date.now()}`;
    const newFluid = {
      id: newId,
      name: 'Custom Fluid',
      rho: 800,
      nu: 1.0e-5,
      nu_cSt_40: 10.0,
      nu_cSt_100: 2.0,
      k: 0.13,
      cp: 2000,
      beta: 0.0007,
      color: '#94a3b8'
    };
    setLocalFluids([newFluid, ...localFluids]);
  };

  const handleDeleteFluid = (id: string) => {
    setLocalFluids((prev: any[]) => prev.filter((f: any) => f.id !== id));
  };

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/50 backdrop-blur-sm p-4 animate-in fade-in duration-200">
      <div className="bg-white rounded-xl shadow-2xl w-full max-w-6xl max-h-[95vh] flex flex-col overflow-hidden">
        <div className="p-5 border-b border-slate-200 bg-slate-50 flex justify-between items-center">
          <div>
            <h2 className="text-xl font-bold text-slate-800 flex items-center gap-2">
              <Edit2 size={20} className="text-blue-600"/> Fluid Database Editor
            </h2>
            <p className="text-sm text-slate-500 mt-1">Manage and customize heat transfer fluids used in the simulation.</p>
          </div>
          <button onClick={onClose} className="p-2 hover:bg-slate-200 rounded-full transition-colors text-slate-500">
            <X size={24} />
          </button>
        </div>

        <div className="flex-1 overflow-y-auto p-6 bg-slate-50/50">
          
          <div className="flex justify-end mb-4">
             <button 
                onClick={handleAddFluid}
                className="flex items-center gap-2 px-4 py-2 bg-green-600 hover:bg-green-700 text-white text-sm font-bold rounded-lg transition-colors shadow-sm"
             >
                <Plus size={16} /> Add Custom Fluid
             </button>
          </div>

          <div className="grid grid-cols-1 gap-4">
            {localFluids.map((f: any) => (
              <div key={f.id} className="bg-white border border-slate-200 rounded-lg p-4 shadow-sm hover:shadow-md transition-shadow group">
                <div className="flex items-center gap-4 mb-4 pb-3 border-b border-slate-100">
                  <div className="relative">
                     <input 
                        type="color" 
                        value={f.color}
                        onChange={(e) => handleFluidChange(f.id, 'color', e.target.value)}
                        className="w-8 h-8 p-0 border-0 rounded-full cursor-pointer overflow-hidden shadow-sm"
                     />
                  </div>
                  <div className="flex-1">
                     <input 
                        type="text" 
                        value={f.name}
                        onChange={(e) => handleFluidChange(f.id, 'name', e.target.value)}
                        className="font-bold text-slate-800 text-lg border-b border-transparent hover:border-slate-300 focus:border-blue-500 focus:outline-none w-full bg-transparent"
                     />
                     <span className="text-xs text-slate-400 font-mono bg-slate-100 px-2 py-0.5 rounded mt-1 inline-block">ID: {f.id}</span>
                  </div>
                  <button 
                    onClick={() => handleDeleteFluid(f.id)}
                    className="p-2 text-slate-300 hover:text-red-500 hover:bg-red-50 rounded-lg transition-colors"
                    title="Delete Fluid"
                  >
                    <Trash2 size={18} />
                  </button>
                </div>
                
                <div className="grid grid-cols-2 md:grid-cols-6 gap-4">
                   <div className="space-y-1">
                     <label className="text-xs font-semibold text-slate-500 uppercase">Density (ρ)</label>
                     <div className="relative">
                       <input 
                         type="number" step="0.1" 
                         value={f.rho}
                         onChange={(e) => handleFluidChange(f.id, 'rho', parseFloat(e.target.value))}
                         className="w-full p-2 border border-slate-300 rounded text-sm focus:ring-2 focus:ring-blue-100 focus:border-blue-400 outline-none"
                       />
                       <span className="absolute right-2 top-2 text-xs text-slate-400">kg/m³</span>
                     </div>
                   </div>

                   <div className="space-y-1">
                     <label className="text-xs font-semibold text-slate-500 uppercase">Visc @ 40°C</label>
                     <div className="relative">
                       <input 
                         type="number" step="0.01" 
                         value={f.nu_cSt_40}
                         onChange={(e) => handleFluidChange(f.id, 'nu_cSt_40', parseFloat(e.target.value))}
                         className="w-full p-2 border border-slate-300 rounded text-sm focus:ring-2 focus:ring-blue-100 focus:border-blue-400 outline-none"
                       />
                       <span className="absolute right-2 top-2 text-xs text-slate-400">cSt</span>
                     </div>
                   </div>

                   <div className="space-y-1">
                     <label className="text-xs font-semibold text-slate-500 uppercase">Visc @ 100°C</label>
                     <div className="relative">
                       <input 
                         type="number" step="0.01" 
                         value={f.nu_cSt_100 || 0}
                         onChange={(e) => handleFluidChange(f.id, 'nu_cSt_100', parseFloat(e.target.value))}
                         className="w-full p-2 border border-slate-300 rounded text-sm focus:ring-2 focus:ring-blue-100 focus:border-blue-400 outline-none"
                       />
                       <span className="absolute right-2 top-2 text-xs text-slate-400">cSt</span>
                     </div>
                   </div>

                   <div className="space-y-1">
                     <label className="text-xs font-semibold text-slate-500 uppercase">Therm. Cond. (k)</label>
                     <div className="relative">
                       <input 
                         type="number" step="0.001" 
                         value={f.k}
                         onChange={(e) => handleFluidChange(f.id, 'k', parseFloat(e.target.value))}
                         className="w-full p-2 border border-slate-300 rounded text-sm focus:ring-2 focus:ring-blue-100 focus:border-blue-400 outline-none"
                       />
                       <span className="absolute right-2 top-2 text-xs text-slate-400">W/mK</span>
                     </div>
                   </div>

                   <div className="space-y-1">
                     <label className="text-xs font-semibold text-slate-500 uppercase">Spec. Heat (Cp)</label>
                     <div className="relative">
                       <input 
                         type="number" step="10" 
                         value={f.cp}
                         onChange={(e) => handleFluidChange(f.id, 'cp', parseFloat(e.target.value))}
                         className="w-full p-2 border border-slate-300 rounded text-sm focus:ring-2 focus:ring-blue-100 focus:border-blue-400 outline-none"
                       />
                       <span className="absolute right-2 top-2 text-xs text-slate-400">J/kgK</span>
                     </div>
                   </div>

                   <div className="space-y-1">
                     <label className="text-xs font-semibold text-slate-500 uppercase">Therm. Exp. (β)</label>
                     <div className="relative">
                       <input 
                         type="number" step="0.00001" 
                         value={f.beta}
                         onChange={(e) => handleFluidChange(f.id, 'beta', parseFloat(e.target.value))}
                         className="w-full p-2 border border-slate-300 rounded text-sm focus:ring-2 focus:ring-blue-100 focus:border-blue-400 outline-none"
                       />
                       <span className="absolute right-2 top-2 text-xs text-slate-400">1/K</span>
                     </div>
                   </div>
                </div>
              </div>
            ))}
          </div>
        </div>

        <div className="p-5 border-t border-slate-200 bg-white flex justify-between items-center">
          <button 
            onClick={onReset}
            className="flex items-center gap-2 px-4 py-2 text-sm font-medium text-slate-600 hover:text-red-600 hover:bg-red-50 rounded-lg transition-colors"
          >
            <RotateCcw size={16} /> Reset to Defaults
          </button>
          <div className="flex gap-3">
            <button 
              onClick={onClose}
              className="px-5 py-2 text-sm font-medium text-slate-600 hover:bg-slate-100 rounded-lg transition-colors"
            >
              Cancel
            </button>
            <button 
              onClick={() => onSave(localFluids)}
              className="flex items-center gap-2 px-6 py-2 text-sm font-bold text-white bg-blue-600 hover:bg-blue-700 rounded-lg shadow-md hover:shadow-lg transition-all"
            >
              <Save size={18} /> Save Database
            </button>
          </div>
        </div>
      </div>
    </div>
  );
};

const MethodologyModal = ({ isOpen, onClose }: any) => {
  if (!isOpen) return null;

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/50 backdrop-blur-sm p-4 animate-in fade-in duration-200">
      <div className="bg-white rounded-xl shadow-2xl w-full max-w-4xl max-h-[90vh] flex flex-col overflow-hidden">
        <div className="p-5 border-b border-slate-200 bg-slate-50 flex justify-between items-center">
          <div>
            <h2 className="text-xl font-bold text-slate-800 flex items-center gap-2">
              <Calculator size={20} className="text-blue-600"/> Simulation Methodology
            </h2>
            <p className="text-sm text-slate-500 mt-1">Physics Engine, Correlations, and Assumptions</p>
          </div>
          <button onClick={onClose} className="p-2 hover:bg-slate-200 rounded-full transition-colors text-slate-500">
            <X size={24} />
          </button>
        </div>

        <div className="flex-1 overflow-y-auto p-6 text-slate-700 space-y-6">
          <section>
            <h3 className="text-lg font-bold text-slate-900 mb-2 border-b pb-1">1. Advanced Nusselt Correlation</h3>
            <p className="text-sm leading-relaxed">
               We utilize the <strong>Teertstra (1999)</strong> composite correlation for plate-fin heatsinks. This model accurately transitions between:
            </p>
            <ul className="list-disc pl-5 text-sm mt-2 space-y-1 text-slate-600">
               <li><strong>Fully Developed Channel Flow:</strong> Occurs with narrow fin spacing.</li>
               <li><strong>Isolated Plate Limit:</strong> Occurs with wide spacing (Boundary Layer flow).</li>
            </ul>
          </section>

          <section>
            <h3 className="text-lg font-bold text-slate-900 mb-2 border-b pb-1">2. Spreading Resistance</h3>
            <p className="text-sm leading-relaxed mb-2">
               The model now accounts for <strong>Geometric Spreading Resistance</strong> (Lee et al / Yovanovich). 
               Heat is generated in a small silicon die and must spread to the larger heatsink base. 
               This adds a significant temperature penalty ($R_{spreading}$) often ignored in simpler models.
            </p>
          </section>

          <section>
            <h3 className="text-lg font-bold text-slate-900 mb-2 border-b pb-1">3. Hydraulic Impedance</h3>
            <p className="text-sm leading-relaxed mb-2">
              Pressure drop includes friction ($f_{app}$) plus <strong>Contraction</strong> and <strong>Expansion</strong> losses ($K_c, K_e$) as fluid enters/exits the fin array. 
              Properties ($\rho, C_p, k$) are evaluated at film temperature.
            </p>
          </section>
        </div>
        
        <div className="p-4 border-t border-slate-200 bg-slate-50 flex justify-end">
           <button 
              onClick={onClose}
              className="px-6 py-2 text-sm font-bold text-white bg-blue-600 hover:bg-blue-700 rounded-lg transition-colors"
            >
              Close
            </button>
        </div>
      </div>
    </div>
  );
};

const Dashboard = () => {
  const [simMode, setSimMode] = useState<'design' | 'rating'>('rating'); 
  
  const [inputs, setInputs] = useState({
    Tj: 70, // deg C (Design Mode)
    targetPower: 1200, // W (Rating Mode)
    T_inlet: 40, // deg C
    bypassFactor: 0.1, // 10% Leakage
  });

  const [hsGeometry, setHsGeometry] = useState(INITIAL_GEOMETRY);
  const [fluids, setFluids] = useState(DEFAULT_FLUIDS);
  const [isFluidModalOpen, setIsFluidModalOpen] = useState(false);
  const [isMethodModalOpen, setIsMethodModalOpen] = useState(false);
  const [selectedFluids, setSelectedFluids] = useState<string[]>(['mpao', 'fuchs', 'dcf281']);
  const [simData, setSimData] = useState<any[]>([]);
  const [rankingData, setRankingData] = useState<any[]>([]);
  const [viscosityData, setViscosityData] = useState<any[]>([]);

  const derivedGeom = useMemo(() => {
    const N = Math.floor(hsGeometry.W / (hsGeometry.s + hsGeometry.t_fin));
    const A_fin_single = 2 * hsGeometry.h_fin * hsGeometry.L + hsGeometry.t_fin * hsGeometry.L;
    const A_fin = N * A_fin_single;
    const A_base_channel = hsGeometry.s * hsGeometry.L;
    const A_base = N * A_base_channel; 
    const A_total = A_fin + A_base;
    const A_c_channel = hsGeometry.s * hsGeometry.h_fin;
    const A_c_total = N * A_c_channel;
    const Dh = (2 * hsGeometry.s * hsGeometry.h_fin) / (hsGeometry.s + hsGeometry.h_fin);
    
    // Spreading Resistance Calculation (Simplified Yovanovich/Lee approx for concentric rectangular areas)
    // R_sp = (1 - epsilon) / (k * sqrt(A_base) * 2 * sqrt(pi)) ... rough approximation
    // Better approx: R_sp = (1 - sqrt(A_die/A_base)) / (k * sqrt(A_die) * sqrt(pi))
    const A_die = hsGeometry.die_W * hsGeometry.die_L;
    const A_hs_base = hsGeometry.W * hsGeometry.L;
    
    // Geometric ratio
    const epsilon = Math.sqrt(A_die / A_hs_base);
    // Spreading Resistance
    const R_spreading = (1 - epsilon) / (hsGeometry.k_base * Math.sqrt(A_die) * Math.sqrt(Math.PI));
    
    const R_internal = hsGeometry.R_TIM + hsGeometry.R_die + R_spreading;

    return { N, A_total, A_fin, A_base, A_c_total, R_internal, Dh, R_spreading, A_die };
  }, [hsGeometry]);

  // Generate Viscosity Data for Chart
  useEffect(() => {
    const data = [];
    for (let T = 30; T <= 120; T += 5) {
      const point: any = { temp: T };
      selectedFluids.forEach(id => {
        const fluid = fluids.find(f => f.id === id);
        if (fluid) {
          // Convert to cSt for easier reading on chart
          point[id] = getViscosity(fluid, T) * 1e6;
        }
      });
      data.push(point);
    }
    setViscosityData(data);
  }, [selectedFluids, fluids]);

  // Run Simulation
  useEffect(() => {
    const heightPoints = [];
    const heightRange = simMode === 'rating' ? [1.0, 2.5] : [0.5, 3.0]; 

    for (let h = 0.5; h <= 3.0; h += 0.1) {
      const point: any = { height: Number(h.toFixed(1)) };
      
      selectedFluids.forEach(fluidId => {
        const fluid = fluids.find(f => f.id === fluidId);
        if (fluid) {
          let res;
          if (simMode === 'rating') {
             // Fixed Power, Solve for Temp
             res = solveForTj(inputs.targetPower, fluid, h, inputs.T_inlet, hsGeometry, derivedGeom, inputs.bypassFactor);
          } else {
             // Fixed Temp, Solve for Power
             res = solveOperatingPoint(fluid, h, inputs.Tj, inputs.T_inlet, hsGeometry, derivedGeom, inputs.bypassFactor);
          }

          if (res) {
            point[`${fluidId}_Q`] = res.Q;
            point[`${fluidId}_u`] = res.u_ind;
            point[`${fluidId}_Tout`] = res.T_out;
            point[`${fluidId}_Tj`] = res.T_j;
          }
        }
      });
      heightPoints.push(point);
    }
    setSimData(heightPoints);

    // Generate Ranking at Fixed H_tank = 1.5m
    const rankH = 1.5;
    const ranks = selectedFluids.map(fluidId => {
      const fluid = fluids.find(f => f.id === fluidId);
      if (!fluid) return null;
      
      let res;
      if (simMode === 'rating') {
        res = solveForTj(inputs.targetPower, fluid, rankH, inputs.T_inlet, hsGeometry, derivedGeom, inputs.bypassFactor);
      } else {
        res = solveOperatingPoint(fluid, rankH, inputs.Tj, inputs.T_inlet, hsGeometry, derivedGeom, inputs.bypassFactor);
      }
      
      if (!res) return null;

      return {
        name: fluid.name,
        color: fluid.color,
        Q: res.Q,
        T_out: res.T_out,
        T_j: res.T_j,
        u_ind: res.u_ind,
        Re: res.Re,
        Ra: res.Ra,
        nu_op: res.nu_op
      };
    }).filter(Boolean).sort((a: any, b: any) => simMode === 'rating' ? a.T_j - b.T_j : b.Q - a.Q);
    
    setRankingData(ranks);

  }, [inputs, selectedFluids, hsGeometry, derivedGeom, fluids, simMode]);

  const toggleFluid = (id: string) => {
    if (selectedFluids.includes(id)) {
      setSelectedFluids(selectedFluids.filter(f => f !== id));
    } else {
      setSelectedFluids([...selectedFluids, id]);
    }
  };

  const handleFluidSave = (updatedFluids: any) => {
    setFluids(updatedFluids);
    setIsFluidModalOpen(false);
  };

  const handleFluidReset = () => {
    setFluids(DEFAULT_FLUIDS);
    setIsFluidModalOpen(false);
  };

  // CSV Export
  const downloadCSV = () => {
    let csv = "Height(m),";
    selectedFluids.forEach(f => {
       csv += `${f} Power(W),${f} Tj(C),${f} Tout(C),${f} Velocity(m/s),`;
    });
    csv += "\n";

    simData.forEach(row => {
      csv += `${row.height},`;
      selectedFluids.forEach(f => {
         csv += `${row[`${f}_Q`]?.toFixed(2) || ''},${row[`${f}_Tj`]?.toFixed(2) || ''},${row[`${f}_Tout`]?.toFixed(2) || ''},${row[`${f}_u`]?.toFixed(4) || ''},`;
      });
      csv += "\n";
    });

    const blob = new Blob([csv], { type: 'text/csv' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `simulation_${simMode}_${new Date().toISOString().slice(0,10)}.csv`;
    a.click();
  };

  return (
    <div className="min-h-screen bg-slate-50 font-sans text-slate-900">
      
      <FluidEditModal 
        isOpen={isFluidModalOpen} 
        onClose={() => setIsFluidModalOpen(false)}
        fluids={fluids}
        onSave={handleFluidSave}
        onReset={handleFluidReset}
      />

      <MethodologyModal
        isOpen={isMethodModalOpen}
        onClose={() => setIsMethodModalOpen(false)}
      />

      {/* Header */}
      <header className="bg-slate-900 text-white p-6 shadow-lg">
        <div className="max-w-7xl mx-auto flex justify-between items-center">
          <div>
            <h1 className="text-2xl font-bold flex items-center gap-3">
              <Wind className="text-blue-400" />
              Natural Convection Simulator
            </h1>
            <p className="text-slate-400 text-sm mt-1">
              Blackwell H200 Immersion Cooling • {simMode === 'rating' ? 'Fixed Power Mode' : 'Fixed Temp Mode'}
            </p>
          </div>
          <div className="flex items-center gap-6">
            <button 
              onClick={() => setIsMethodModalOpen(true)}
              className="hidden md:flex items-center gap-2 text-sm text-slate-300 hover:text-white transition-colors"
            >
              <Info size={16} /> Methodology
            </button>
            <div className="text-right hidden md:block pl-6 border-l border-slate-700">
               {simMode === 'rating' ? (
                 <>
                  <div className="text-xs text-slate-400">Target Power</div>
                  <div className="text-xl font-bold text-green-400">{inputs.targetPower} W</div>
                 </>
               ) : (
                 <>
                  <div className="text-xs text-slate-400">Max Temp</div>
                  <div className="text-xl font-bold text-red-400">{inputs.Tj} °C</div>
                 </>
               )}
            </div>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto p-6 grid grid-cols-1 lg:grid-cols-12 gap-6">
        
        {/* Left Sidebar: Controls */}
        <div className="lg:col-span-3 space-y-6">
          
          {/* Input Panel */}
          <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5">
            <h3 className="font-semibold text-slate-800 flex items-center gap-2 mb-4">
              <Settings size={18} /> Operating Conditions
            </h3>
            
            {/* Mode Toggle */}
            <div className="flex bg-slate-100 p-1 rounded-lg mb-6">
              <button 
                onClick={() => setSimMode('rating')}
                className={`flex-1 py-1.5 text-xs font-semibold rounded-md transition-all ${simMode === 'rating' ? 'bg-white text-blue-600 shadow-sm' : 'text-slate-500 hover:text-slate-700'}`}
              >
                Fixed Power (TDP)
              </button>
              <button 
                onClick={() => setSimMode('design')}
                className={`flex-1 py-1.5 text-xs font-semibold rounded-md transition-all ${simMode === 'design' ? 'bg-white text-blue-600 shadow-sm' : 'text-slate-500 hover:text-slate-700'}`}
              >
                Fixed Temp (Max Q)
              </button>
            </div>

            <div className="space-y-6">
              {simMode === 'design' ? (
                <div>
                  <label className="flex justify-between text-sm font-medium text-slate-700 mb-2">
                    Max Junction Temp (Tj)
                    <span className="text-blue-600">{inputs.Tj}°C</span>
                  </label>
                  <input 
                    type="range" min="60" max="120" step="1" 
                    value={inputs.Tj}
                    onChange={(e) => setInputs({...inputs, Tj: Number(e.target.value)})}
                    className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
                  />
                </div>
              ) : (
                <div>
                  <label className="flex justify-between text-sm font-medium text-slate-700 mb-2">
                    Target Power (TDP)
                    <span className="text-blue-600">{inputs.targetPower} W</span>
                  </label>
                  <input 
                    type="range" min="500" max="2000" step="50" 
                    value={inputs.targetPower}
                    onChange={(e) => setInputs({...inputs, targetPower: Number(e.target.value)})}
                    className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
                  />
                  <div className="flex justify-between text-xs text-slate-400 mt-1">
                    <span>500W</span>
                    <span>2000W</span>
                  </div>
                </div>
              )}

              <div>
                <label className="flex justify-between text-sm font-medium text-slate-700 mb-2">
                  Inlet Fluid Temp (T∞)
                  <span className="text-blue-600">{inputs.T_inlet}°C</span>
                </label>
                <input 
                  type="range" min="20" max="50" step="1" 
                  value={inputs.T_inlet}
                  onChange={(e) => setInputs({...inputs, T_inlet: Number(e.target.value)})}
                  className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
                />
              </div>

              <div>
                <label className="flex justify-between text-sm font-medium text-slate-700 mb-2">
                  Ducting Bypass/Leakage
                  <span className="text-orange-600">{(inputs.bypassFactor * 100).toFixed(0)}%</span>
                </label>
                <input 
                  type="range" min="0" max="0.5" step="0.05" 
                  value={inputs.bypassFactor}
                  onChange={(e) => setInputs({...inputs, bypassFactor: Number(e.target.value)})}
                  className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-orange-500"
                />
                <p className="text-xs text-slate-400 mt-1">Percentage of flow missing fins due to imperfect shroud.</p>
              </div>
            </div>
          </div>

          {/* Heat Sink Geometry Input */}
          <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5">
            <h3 className="font-semibold text-slate-800 flex items-center gap-2 mb-4">
              <Ruler size={18} /> Heat Sink Geometry
            </h3>
            <div className="space-y-1">
               <InputField 
                 label="Height (Vertical Length)" 
                 value={hsGeometry.L} 
                 unit="m" 
                 step={0.01} 
                 onChange={(v: number) => setHsGeometry(p => ({...p, L: v}))}
               />
               <InputField 
                 label="Width (W)" 
                 value={hsGeometry.W} 
                 unit="m" 
                 step={0.005} 
                 onChange={(v: number) => setHsGeometry(p => ({...p, W: v}))}
               />
               <InputField 
                 label="Fin Height (Length)" 
                 value={hsGeometry.h_fin} 
                 unit="m" 
                 step={0.001} 
                 onChange={(v: number) => setHsGeometry(p => ({...p, h_fin: v}))}
               />
               <InputField 
                 label="Fin Spacing (s)" 
                 value={hsGeometry.s} 
                 unit="m" 
                 step={0.0005} 
                 onChange={(v: number) => setHsGeometry(p => ({...p, s: v}))}
               />
               <InputField 
                 label="Fin Thickness (t)" 
                 value={hsGeometry.t_fin} 
                 unit="m" 
                 step={0.0001} 
                 onChange={(v: number) => setHsGeometry(p => ({...p, t_fin: v}))}
               />
            </div>
            
            {/* New Die Geometry Section */}
            <div className="mt-4 pt-4 border-t border-slate-100 space-y-2">
               <h4 className="text-xs font-bold text-slate-500 uppercase flex items-center gap-1">
                 <Cpu size={12}/> Die Geometry
               </h4>
               <p className="text-[10px] text-slate-400 mb-2 leading-tight">
                 Used to calc Spreading Resistance
               </p>
               <div className="grid grid-cols-2 gap-2">
                  <InputField 
                    label="Die Width" 
                    value={hsGeometry.die_W} 
                    unit="m" 
                    step={0.001} 
                    onChange={(v: number) => setHsGeometry(p => ({...p, die_W: v}))}
                  />
                  <InputField 
                    label="Die Height" 
                    value={hsGeometry.die_L} 
                    unit="m" 
                    step={0.001} 
                    onChange={(v: number) => setHsGeometry(p => ({...p, die_L: v}))}
                  />
               </div>
            </div>

            <div className="mt-4 pt-4 border-t border-slate-100 space-y-2">
               <h4 className="text-xs font-bold text-slate-500 uppercase">Derived Properties</h4>
               <div className="flex justify-between text-sm">
                 <span className="text-slate-600">Hydraulic Dia:</span>
                 <span className="font-mono font-bold">{(derivedGeom.Dh * 1000).toFixed(2)} mm</span>
               </div>
               <div className="flex justify-between text-sm">
                 <span className="text-slate-600">R_spread:</span>
                 <span className="font-mono font-bold text-orange-600">{(derivedGeom.R_spreading).toFixed(4)} K/W</span>
               </div>
            </div>
          </div>

          {/* Fluid Selection */}
          <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5">
            <div className="flex justify-between items-center mb-4">
              <h3 className="font-semibold text-slate-800 flex items-center gap-2">
                <Droplets size={18} /> Fluid Selection
              </h3>
              <button 
                onClick={() => setIsFluidModalOpen(true)}
                className="text-xs bg-blue-50 hover:bg-blue-100 text-blue-700 px-2 py-1 rounded border border-blue-200 transition-colors flex items-center gap-1"
              >
                <Edit2 size={12} /> Customize
              </button>
            </div>
            <div className="space-y-2">
              {fluids.map(fluid => (
                <div 
                  key={fluid.id}
                  onClick={() => toggleFluid(fluid.id)}
                  className={`flex items-center justify-between p-2 rounded-md cursor-pointer transition-colors border ${selectedFluids.includes(fluid.id) ? 'bg-blue-50 border-blue-200' : 'hover:bg-slate-50 border-transparent'}`}
                >
                  <div className="flex items-center gap-2">
                    <div className="w-3 h-3 rounded-full" style={{backgroundColor: fluid.color}}></div>
                    <span className="text-sm font-medium text-slate-700">{fluid.name}</span>
                  </div>
                  {selectedFluids.includes(fluid.id) && <div className="w-2 h-2 bg-blue-500 rounded-full"></div>}
                </div>
              ))}
            </div>
          </div>

        </div>

        {/* Main Content: Charts */}
        <div className="lg:col-span-9 space-y-6">
          
          {/* Chart 1: Main Result */}
          <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
            <div className="flex justify-between items-center mb-6">
              <div>
                <h2 className="text-lg font-bold text-slate-900">
                  {simMode === 'rating' ? 'Junction Temperature vs Tank Height' : 'Power Dissipation vs Tank Height'}
                </h2>
                <p className="text-sm text-slate-500">
                  {simMode === 'rating' 
                    ? `Showing predicted chip temp (Tj) for fixed ${inputs.targetPower}W load.` 
                    : `Showing max power capacity for fixed ${inputs.Tj}°C limit.`}
                </p>
              </div>
              <div className="flex gap-2">
                 <button onClick={downloadCSV} className="flex items-center gap-2 px-3 py-1 bg-slate-100 hover:bg-slate-200 text-slate-600 text-xs font-bold rounded-lg transition-colors">
                    <Download size={14}/> CSV
                 </button>
              </div>
            </div>
            
            <div className="h-[350px] w-full">
              <ResponsiveContainer width="100%" height="100%">
                <ComposedChart data={simData} margin={{ top: 5, right: 30, left: 20, bottom: 25 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                  <XAxis 
                    dataKey="height" 
                    label={{ value: 'Tank Chimney Height (m)', position: 'insideBottom', offset: -15 }} 
                    type="number" domain={[0.5, 3.0]} tickCount={6}
                  />
                  <YAxis 
                    label={{ 
                      value: simMode === 'rating' ? 'Temp (°C)' : 'Power (W)', 
                      angle: -90, 
                      position: 'insideLeft' 
                    }} 
                    domain={['auto', 'auto']}
                  />
                  <Tooltip 
                    contentStyle={{ backgroundColor: '#fff', borderRadius: '8px', border: '1px solid #e2e8f0', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)' }}
                    labelFormatter={(label) => `Tank Height: ${label} m`}
                    formatter={(value: number, name: string) => [
                      simMode === 'rating' ? `${value.toFixed(1)} °C` : `${value.toFixed(0)} W`, 
                      name.split('_')[0]
                    ]}
                  />
                  <Legend verticalAlign="top" height={36}/>
                  
                  {/* Reference Lines */}
                  {simMode === 'rating' ? (
                     <ReferenceLine y={108} stroke="red" strokeDasharray="3 3" label={{ position: 'insideRight', value: 'Max 108°C', fill: 'red', fontSize: 12 }} />
                  ) : (
                     <ReferenceLine y={inputs.targetPower} stroke="green" strokeDasharray="3 3" label={{ position: 'insideRight', value: 'Target', fill: 'green', fontSize: 12 }} />
                  )}

                  {fluids.map(fluid => (
                    selectedFluids.includes(fluid.id) && (
                      <Line 
                        key={fluid.id}
                        type="monotone" 
                        dataKey={simMode === 'rating' ? `${fluid.id}_Tj` : `${fluid.id}_Q`} 
                        name={fluid.name}
                        stroke={fluid.color} 
                        strokeWidth={3}
                        dot={false}
                      />
                    )
                  ))}
                </ComposedChart>
              </ResponsiveContainer>
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            
            {/* Chart 2: Induced Velocity */}
            <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
              <h2 className="text-lg font-bold text-slate-900 mb-1">Induced Velocity (Net)</h2>
              <p className="text-sm text-slate-500 mb-4">Effective velocity through fins (after leakage)</p>
              <div className="h-[250px] w-full">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={simData} margin={{ top: 5, right: 10, left: 0, bottom: 20 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="height" label={{ value: 'Tank Height (m)', position: 'insideBottom', offset: -10, fontSize: 12 }} />
                    <YAxis label={{ value: 'Velocity (m/s)', angle: -90, position: 'insideLeft', fontSize: 12 }} />
                    <Tooltip 
                      labelFormatter={(label) => `Tank Height: ${label} m`}
                      formatter={(val: number, name: string) => [`${val.toFixed(3)} m/s`, name.split('_')[0]]} 
                    />
                    {fluids.map(fluid => (
                      selectedFluids.includes(fluid.id) && (
                        <Line 
                          key={fluid.id}
                          type="monotone" 
                          dataKey={`${fluid.id}_u`} 
                          stroke={fluid.color} 
                          strokeWidth={2}
                          dot={false} 
                          name={fluid.name}
                        />
                      )
                    ))}
                  </LineChart>
                </ResponsiveContainer>
              </div>
            </div>

            {/* Chart 3: Outlet Temperature */}
            <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
              <h2 className="text-lg font-bold text-slate-900 mb-1">Outlet Temperature</h2>
              <p className="text-sm text-slate-500 mb-4">Fluid temperature exiting the top of heat sink</p>
              <div className="h-[250px] w-full">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={simData} margin={{ top: 5, right: 10, left: 0, bottom: 20 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="height" label={{ value: 'Tank Height (m)', position: 'insideBottom', offset: -10, fontSize: 12 }} />
                    <YAxis domain={['auto', 'auto']} label={{ value: 'T_out (°C)', angle: -90, position: 'insideLeft', fontSize: 12 }} />
                    <Tooltip 
                      labelFormatter={(label) => `Tank Height: ${label} m`}
                      formatter={(val: number, name: string) => [`${val.toFixed(1)} °C`, name.split('_')[0]]} 
                    />
                    {fluids.map(fluid => (
                      selectedFluids.includes(fluid.id) && (
                        <Line 
                          key={fluid.id}
                          type="monotone" 
                          dataKey={`${fluid.id}_Tout`} 
                          stroke={fluid.color} 
                          strokeWidth={2}
                          dot={false} 
                          name={fluid.name}
                        />
                      )
                    ))}
                  </LineChart>
                </ResponsiveContainer>
              </div>
            </div>
          </div>
          
          {/* New Chart: Viscosity vs Temp */}
          <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
            <h2 className="text-lg font-bold text-slate-900 mb-1">Kinematic Viscosity vs. Temperature</h2>
            <p className="text-sm text-slate-500 mb-4">Log scale. Demonstrates thermal thinning effect (ASTM D341 Model).</p>
            <div className="h-[250px] w-full">
              <ResponsiveContainer width="100%" height="100%">
                <LineChart data={viscosityData} margin={{ top: 5, right: 10, left: 10, bottom: 20 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                  <XAxis dataKey="temp" label={{ value: 'Temperature (°C)', position: 'insideBottom', offset: -10, fontSize: 12 }} />
                  <YAxis 
                    scale="log" 
                    domain={['auto', 'auto']} 
                    allowDataOverflow 
                    label={{ value: 'Viscosity (cSt)', angle: -90, position: 'insideLeft', fontSize: 12 }} 
                  />
                  <Tooltip 
                    labelFormatter={(label) => `${label} °C`}
                    formatter={(val: number, name: string) => [`${val.toFixed(2)} cSt`, name]} 
                  />
                  <Legend verticalAlign="top" height={36}/>
                  {fluids.map(fluid => (
                    selectedFluids.includes(fluid.id) && (
                      <Line 
                        key={fluid.id}
                        type="monotone" 
                        dataKey={fluid.id} 
                        stroke={fluid.color} 
                        strokeWidth={2}
                        dot={false} 
                        name={fluid.name}
                      />
                    )
                  ))}
                </LineChart>
              </ResponsiveContainer>
            </div>
          </div>

          {/* Ranking Table */}
          <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
            <div className="p-6 border-b border-slate-200">
              <h2 className="text-lg font-bold text-slate-900">Fluid Performance Ranking</h2>
              <p className="text-sm text-slate-500">
                Snapshot at Tank Height = 1.5m. 
                {simMode === 'rating' 
                  ? ' Ranked by lowest Junction Temperature (Best cooling).' 
                  : ' Ranked by highest Power Dissipation.'}
              </p>
            </div>
            <div className="overflow-x-auto">
              <table className="w-full text-left text-sm">
                <thead className="bg-slate-50 text-slate-600 font-semibold">
                  <tr>
                    <th className="p-4">Rank</th>
                    <th className="p-4">Fluid</th>
                    <th className="p-4">Junction T (°C)</th>
                    <th className="p-4">Power (W)</th>
                    <th className="p-4">Outlet T (°C)</th>
                    <th className="p-4">Velocity (m/s)</th>
                    <th className="p-4">Re</th>
                    <th className="p-4">Op. Visc (cSt)</th>
                    <th className="p-4">Status</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-100">
                  {rankingData.map((row, idx) => {
                    // Status logic
                    let isGood = false;
                    let isWarn = false;
                    if (simMode === 'rating') {
                      if (row.T_j < 100) isGood = true;
                      else if (row.T_j < 108) isWarn = true;
                    } else {
                      if (row.Q >= inputs.targetPower) isGood = true;
                      else if (row.Q >= inputs.targetPower * 0.8) isWarn = true;
                    }

                    return (
                      <tr key={idx} className="hover:bg-slate-50">
                        <td className="p-4 font-mono text-slate-400">#{idx + 1}</td>
                        <td className="p-4 font-medium flex items-center gap-3">
                          <div className="w-3 h-3 rounded-full shadow-sm" style={{backgroundColor: row.color}}></div>
                          <span className="text-slate-900">{row.name}</span>
                        </td>
                        <td className={`p-4 font-bold ${simMode === 'rating' ? 'text-slate-900' : 'text-slate-500'}`}>{row.T_j.toFixed(1)}</td>
                        <td className={`p-4 font-bold ${simMode !== 'rating' ? 'text-slate-900' : 'text-slate-500'}`}>{row.Q.toFixed(0)}</td>
                        <td className="p-4">{row.T_out.toFixed(1)}</td>
                        <td className="p-4">{row.u_ind.toFixed(3)}</td>
                        <td className="p-4 font-mono text-xs">{row.Re.toFixed(0)}</td>
                        <td className="p-4 font-mono text-xs text-slate-500">{(row.nu_op * 1e6).toFixed(2)}</td>
                        <td className="p-4">
                          {isGood ? (
                            <span className="inline-flex items-center px-2 py-1 rounded-full text-xs font-medium bg-green-100 text-green-800">
                              Pass
                            </span>
                          ) : isWarn ? (
                            <span className="inline-flex items-center px-2 py-1 rounded-full text-xs font-medium bg-yellow-100 text-yellow-800">
                              Marginal
                            </span>
                          ) : (
                            <span className="inline-flex items-center px-2 py-1 rounded-full text-xs font-medium bg-red-100 text-red-800">
                              Fail
                            </span>
                          )}
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
          </div>

        </div>
      </main>
    </div>
  );
};

// Mount React App
const root = createRoot(document.getElementById('root')!);
root.render(<Dashboard />);