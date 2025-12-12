import React, { useState, useEffect, useMemo } from 'react';
import { createRoot } from 'react-dom/client';
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend,
  ResponsiveContainer, ReferenceLine, ComposedChart
} from 'recharts';
import {
  Settings, Wind, Info, Droplets, Ruler, Calculator, Edit2, X,
  Save, RotateCcw, Download, Trash2, Plus, Cpu
} from 'lucide-react';

/* ==============================================================================================
   SECTION 1: TYPES & INTERFACES
   ============================================================================================== */

interface Fluid {
  id: string;
  name: string;
  color: string;
  rho: number;        // Density [kg/m^3]
  k: number;          // Thermal Conductivity [W/m-K]
  cp: number;         // Specific Heat Capacity [J/kg-K]
  beta: number;       // Volumetric Thermal Expansion Coefficient [1/K]
  nu: number;         // Kinematic Viscosity at Reference Temp [m^2/s]
  nu_cSt_40?: number; // Viscosity at 40°C for ASTM D341 [cSt]
  nu_cSt_100?: number;// Viscosity at 100°C for ASTM D341 [cSt]
}

interface Geometry {
  L: number;      // Heatsink Height (Flow Length) [m]
  W: number;      // Heatsink Width [m]
  s: number;      // Fin Spacing [m]
  t_fin: number;  // Fin Thickness [m]
  h_fin: number;  // Fin Extrusion Length (Depth) [m]
  die_W: number;  // Die Width [m]
  die_L: number;  // Die Height [m]
  R_TIM: number;  // Thermal Interface Material Resistance [K/W]
  R_die: number;  // Junction-to-Case Resistance [K/W]
  k_base: number; // Thermal Conductivity of Base Material [W/m-K]
}

interface SimResult {
  Q: number;      // Dissipated Power [W]
  u_ind: number;  // Induced Velocity [m/s]
  T_out: number;  // Outlet Temperature [°C]
  T_wall: number; // Base (Wall) Temperature [°C]
  T_j: number;    // Junction Temperature [°C]
  h_conv: number; // Convection Coefficient [W/m^2-K]
  Re: number;     // Reynolds Number
  Ra: number;     // Rayleigh Number
}

/* ==============================================================================================
   SECTION 2: CONSTANTS & DEFAULTS
   ============================================================================================== */

const CONSTANTS = {
  G: 9.81,          // Gravity [m/s^2]
  T_REF: 273.15,    // Kelvin Offset
  P_ATM: 101325,    // Standard Pressure [Pa]
};

// Default Heatsink Geometry (High-performance Plate Fin Stack)
const INITIAL_GEOMETRY: Geometry = {
  L: 0.15, 
  W: 0.08, 
  s: 0.003, 
  t_fin: 0.0008, 
  h_fin: 0.025,
  die_W: 0.02, 
  die_L: 0.03,
  R_TIM: 0.02, 
  R_die: 0.01, 
  k_base: 400
};

// Default Fluid Database (Single-Phase Dielectric Fluids)
const DEFAULT_FLUIDS: Fluid[] = [
  { 
    id: 'dcf281', name: 'DCF-281 (Novvi)', color: '#64748b',
    rho: 793.36, k: 0.14059, cp: 2110, beta: 0.00092,
    nu: 10.41e-6, nu_cSt_40: 10.41, nu_cSt_100: 2.81 
  },
  { 
    id: 'shell', name: 'Shell XHVI3', color: '#f59e0b',
    rho: 790.2, k: 0.13702, cp: 2070, beta: 0.00064,
    nu: 9.99e-6, nu_cSt_40: 9.99, nu_cSt_100: 2.72 
  },
  { 
    id: 'castrol', name: 'Castrol DC15', color: '#10b981',
    rho: 819, k: 0.134, cp: 2200, beta: 0.00090,
    nu: 7.50e-6, nu_cSt_40: 7.50, nu_cSt_100: 2.20 
  },
  { 
    id: 'valvoline', name: 'Valvoline HPC', color: '#ef4444',
    rho: 811.5, k: 0.1304, cp: 2000, beta: 0.00080,
    nu: 8.01e-6, nu_cSt_40: 8.01, nu_cSt_100: 2.40 
  },
  { 
    id: 'fuchs', name: 'Fuchs Renolin', color: '#8b5cf6',
    rho: 826, k: 0.134, cp: 2200, beta: 0.00065,
    nu: 4.96e-6, nu_cSt_40: 4.96, nu_cSt_100: 1.70 
  },
  { 
    id: 'mpao', name: 'Novel MPAO', color: '#3b82f6',
    rho: 794.8, k: 0.15203, cp: 2270, beta: 0.00080,
    nu: 7.994e-6, nu_cSt_40: 7.994, nu_cSt_100: 1.885 
  }, 
];

/* ==============================================================================================
   SECTION 3: PHYSICS ENGINE
   ============================================================================================== */

/**
 * ASTM D341 Viscosity-Temperature Model
 * Calculates kinematic viscosity based on two reference points (40°C and 100°C).
 * Essential for oils where viscosity drops exponentially with temperature.
 */
const getViscosity = (fluid: Fluid, T_C: number): number => {
  if (!fluid.nu_cSt_100 || !fluid.nu_cSt_40) return fluid.nu;

  const T_K = T_C + CONSTANTS.T_REF;
  const T1 = 40 + CONSTANTS.T_REF;
  const T2 = 100 + CONSTANTS.T_REF;

  // Double-Log Equation: log10(log10(v + 0.7)) = A - B * log10(T)
  const v1 = fluid.nu_cSt_40;
  const v2 = fluid.nu_cSt_100;

  const Z1 = Math.log10(Math.log10(v1 + 0.7));
  const Z2 = Math.log10(Math.log10(v2 + 0.7));

  const B = (Z1 - Z2) / (Math.log10(T2) - Math.log10(T1));
  const A = Z1 + B * Math.log10(T1);

  const Z_target = A - B * Math.log10(T_K);
  const v_cSt = Math.pow(10, Math.pow(10, Z_target)) - 0.7;

  return Math.max(0.1, v_cSt) * 1e-6; // Convert cSt to m²/s
};

/**
 * Teertstra (1999) Composite Correlation
 * Blends Channel Flow (fully developed) and Boundary Layer (developing) regimes.
 * Nu_composite = [ (Nu_fd)^-n + (Nu_dev)^-n ] ^ -1/n
 */
const calculateTeertstraNu = (Re_channel: number, Pr: number): number => {
  // Term 1: Fully Developed Channel Flow (Laminar limit between plates)
  const term1 = Math.pow(Re_channel * Pr * 0.5, -3);
  
  // Term 2: Developing Boundary Layer Flow (Flat plate limit)
  const term2 = Math.pow(0.664 * Math.sqrt(Re_channel) * Math.pow(Pr, 0.33), -3);
  
  // Composite blending
  return Math.pow(term1 + term2, -1/3);
};

/**
 * Solves the Hydraulic & Thermal Operating Point
 * Iteratively balances Buoyancy Pressure vs. Friction/Form Losses.
 */
const solveOperatingPoint = (
  fluid: Fluid, 
  H_tank: number, 
  Tj: number, 
  T_inlet: number,
  geom: Geometry,
  derived: { N: number, A_total: number, A_fin: number, A_c_total: number, R_internal: number, Dh: number, A_base: number },
  bypassFactor: number = 0 
): SimResult => {
  
  // -- 1. Initialization --
  let T_wall = Tj - 5; // Guess base temp
  let T_out = T_inlet + 5; // Guess outlet temp
  let u_induced = 0.05; // Guess velocity
  
  // Hydraulic areas
  const A_frontal = geom.W * geom.h_fin; 
  const sigma = derived.A_c_total / A_frontal; // Contraction ratio
  const Kc = 0.42 * (1 - Math.pow(sigma, 2)); // Contraction loss coeff
  const Ke = Math.pow(1 - sigma, 2);          // Expansion loss coeff
  const K_form = Kc + Ke;

  let converged = false;
  let iterations = 0;
  
  // Result holders
  let Q_final = 0;
  let h_final = 0;
  let Re_final = 0;
  let Ra_final = 0;
  let nu_op = fluid.nu;

  while (!converged && iterations < 100) {
    const T_mean = (T_inlet + T_out) / 2;
    
    // -- 2. Fluid Properties at Film Temp --
    const nu_bulk = getViscosity(fluid, T_mean);
    const nu_wall = getViscosity(fluid, T_wall);
    nu_op = nu_bulk;

    const rho_bulk = fluid.rho * (1 - fluid.beta * (T_mean - 40));
    const rho_inlet = fluid.rho * (1 - fluid.beta * (T_inlet - 40));
    const cp_fluid = fluid.cp * (1 + 0.0015 * (T_mean - 40)); // Linear Cp increase approx

    // -- 3. Buoyancy (Driving Pressure) --
    // dP_buoy = rho * g * beta * dT * Height
    const dT_core = Math.max(0.1, T_mean - T_inlet);
    const dT_chimney = Math.max(0.1, T_out - T_inlet);
    const H_chimney = Math.max(0, H_tank - geom.L);
    
    const dP_buoy = CONSTANTS.G * fluid.beta * rho_inlet * (dT_core * geom.L + dT_chimney * H_chimney);

    // -- 4. Flow Resistance (Resisting Pressure) --
    const u_guess = Math.max(0.001, u_induced);
    const Re_dh = (u_guess * derived.Dh) / nu_bulk;
    
    // Friction factor (Laminar parallel plates f*Re = 24)
    const f_fd = 24 / Re_dh;
    const f_corr = f_fd * Math.pow(nu_wall / nu_bulk, 0.58); // Sieder-Tate viscosity correction
    
    const K_friction = f_corr * (geom.L / derived.Dh);
    const K_total = K_friction + K_form + 1.0; // +1.0 for exit velocity head loss

    // -- 5. Velocity Solution (Bernoulli) --
    let u_new = Math.sqrt( (2 * dP_buoy) / (rho_bulk * K_total) );
    u_new = u_new * (1 - bypassFactor); // Apply leakage penalty
    
    // Relaxation
    u_induced = u_induced * 0.6 + u_new * 0.4;

    // -- 6. Heat Transfer (Teertstra) --
    const Re = (u_induced * derived.Dh) / nu_bulk;
    const Pr = (nu_bulk * rho_bulk * cp_fluid) / fluid.k;
    
    const Nu = calculateTeertstraNu(Re, Pr);
    const h_bare = (Nu * fluid.k) / derived.Dh;
    
    // Fin Efficiency (Adiabatic tip)
    const m = Math.sqrt((2 * h_bare) / (geom.k_base * geom.t_fin));
    const eta_fin = Math.tanh(m * geom.h_fin) / (m * geom.h_fin);
    
    // Overall Surface Efficiency
    const eta_o = 1 - (derived.A_fin / derived.A_total) * (1 - eta_fin);
    const h_eff = h_bare * eta_o;
    
    // -- 7. Energy Balance --
    const dT_log_mean = Math.max(0.1, T_wall - T_mean); // Simplified LMTD
    const Q_conv = h_eff * derived.A_total * dT_log_mean;
    
    const m_dot = rho_bulk * u_induced * derived.A_c_total;
    const T_out_calc = T_inlet + (Q_conv / (m_dot * cp_fluid));
    
    // Internal Resistance (Junction to Wall)
    const T_wall_calc = Tj - Q_conv * derived.R_internal;

    // -- 8. Convergence Check --
    if (Math.abs(T_out_calc - T_out) < 0.05 && Math.abs(T_wall_calc - T_wall) < 0.05) {
      converged = true;
      Q_final = Q_conv;
      T_out = T_out_calc;
      T_wall = T_wall_calc;
      h_final = h_eff;
      Re_final = Re;
      Ra_final = (CONSTANTS.G * fluid.beta * (T_wall - T_inlet) * Math.pow(geom.s, 3)) / (nu_bulk**2) * Pr;
    } else {
      T_out = T_out * 0.8 + T_out_calc * 0.2;
      T_wall = T_wall * 0.8 + T_wall_calc * 0.2;
      
      // Safety bounds
      if (T_wall > Tj) T_wall = Tj;
      if (T_wall < T_inlet) T_wall = T_inlet + 0.1;
    }
    iterations++;
  }

  return { Q: Q_final, u_ind: u_induced, T_out, T_wall, T_j: Tj, h_conv: h_final, Re: Re_final, Ra: Ra_final };
};

/**
 * Inverse Solver: Finds Junction Temp (Tj) required to dissipate Target Power (Q).
 * Uses Binary Search wrapper around the main solver.
 */
const solveForTj = (
  targetQ: number,
  fluid: Fluid, 
  H_tank: number, 
  T_inlet: number,
  geom: Geometry,
  derived: any,
  bypassFactor: number
): SimResult | null => {
  let low = T_inlet + 0.1;
  let high = 200; // Max safety limit
  let bestResult = null;
  let minError = Infinity;

  for(let i=0; i<20; i++) {
     const midTj = (low + high) / 2;
     const res = solveOperatingPoint(fluid, H_tank, midTj, T_inlet, geom, derived, bypassFactor);
     
     const diff = res.Q - targetQ;
     if (Math.abs(diff) < minError) {
       minError = Math.abs(diff);
       bestResult = res;
     }

     if (Math.abs(diff) < 0.5) break; 
     
     if (res.Q > targetQ) high = midTj; // Too much cooling, lower temp
     else low = midTj; // Not enough cooling, raise temp
  }
  return bestResult;
};

/* ==============================================================================================
   SECTION 4: UI COMPONENTS
   ============================================================================================== */

const InputField = ({ label, value, unit, step = 0.001, min = 0, onChange }: any) => (
  <div className="mb-3">
    <label className="flex justify-between text-xs font-semibold text-slate-600 mb-1 uppercase tracking-wider">
      {label}
    </label>
    <div className="flex items-center relative">
      <input
        type="number" step={step} min={min} value={value}
        onChange={(e) => onChange(parseFloat(e.target.value) || 0)}
        className="w-full px-3 py-2 bg-white border border-slate-300 rounded text-sm text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 font-mono"
      />
      <span className="absolute right-3 text-xs text-slate-400 pointer-events-none">{unit}</span>
    </div>
  </div>
);

const FluidEditModal = ({ isOpen, onClose, fluids, onSave, onReset }: any) => {
  const [localFluids, setLocalFluids] = useState<Fluid[]>(fluids);

  useEffect(() => {
    if (isOpen) setLocalFluids(JSON.parse(JSON.stringify(fluids)));
  }, [fluids, isOpen]);

  if (!isOpen) return null;

  const handleFluidChange = (id: string, field: keyof Fluid, val: any) => {
    setLocalFluids(prev => prev.map(f => f.id === id ? { ...f, [field]: val } : f));
  };

  const handleAddFluid = () => {
    setLocalFluids([{
      id: `custom_${Date.now()}`, name: 'New Fluid', color: '#94a3b8',
      rho: 800, k: 0.13, nu: 1.0e-5, nu_cSt_40: 10.0, nu_cSt_100: 2.0, cp: 2000, beta: 0.0007
    }, ...localFluids]);
  };

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/50 backdrop-blur-sm p-4 animate-in fade-in duration-200">
      <div className="bg-white rounded-xl shadow-2xl w-full max-w-6xl max-h-[95vh] flex flex-col overflow-hidden">
        <div className="p-5 border-b border-slate-200 bg-slate-50 flex justify-between items-center">
          <h2 className="text-xl font-bold text-slate-800 flex items-center gap-2">
            <Edit2 size={20} className="text-blue-600"/> Fluid Database Editor
          </h2>
          <button onClick={onClose} className="p-2 hover:bg-slate-200 rounded-full"><X size={24} /></button>
        </div>
        <div className="flex-1 overflow-y-auto p-6 bg-slate-50/50">
          <div className="flex justify-end mb-4">
             <button onClick={handleAddFluid} className="flex items-center gap-2 px-4 py-2 bg-green-600 hover:bg-green-700 text-white text-sm font-bold rounded-lg shadow-sm">
                <Plus size={16} /> Add Custom Fluid
             </button>
          </div>
          <div className="grid grid-cols-1 gap-4">
            {localFluids.map((f) => (
              <div key={f.id} className="bg-white border border-slate-200 rounded-lg p-4 shadow-sm">
                <div className="flex items-center gap-4 mb-4 pb-3 border-b border-slate-100">
                  <input type="color" value={f.color} onChange={(e) => handleFluidChange(f.id, 'color', e.target.value)} className="w-8 h-8 rounded-full cursor-pointer"/>
                  <div className="flex-1">
                     <input type="text" value={f.name} onChange={(e) => handleFluidChange(f.id, 'name', e.target.value)} className="font-bold text-slate-800 text-lg w-full bg-transparent border-b border-transparent hover:border-slate-300 focus:outline-none"/>
                  </div>
                  <button onClick={() => setLocalFluids(prev => prev.filter(lf => lf.id !== f.id))} className="text-slate-300 hover:text-red-500"><Trash2 size={18} /></button>
                </div>
                <div className="grid grid-cols-2 md:grid-cols-6 gap-4">
                   {[
                     { l: 'Density', k: 'rho', u: 'kg/m³', s: 0.1 },
                     { l: 'Visc @ 40C', k: 'nu_cSt_40', u: 'cSt', s: 0.01 },
                     { l: 'Visc @ 100C', k: 'nu_cSt_100', u: 'cSt', s: 0.01 },
                     { l: 'Therm Cond', k: 'k', u: 'W/mK', s: 0.001 },
                     { l: 'Spec Heat', k: 'cp', u: 'J/kgK', s: 10 },
                     { l: 'Expansion', k: 'beta', u: '1/K', s: 0.00001 }
                   ].map((field) => (
                     <div key={field.k} className="space-y-1">
                       <label className="text-xs font-semibold text-slate-500 uppercase">{field.l}</label>
                       <div className="relative">
                         {/* @ts-ignore */}
                         <input type="number" step={field.s} value={f[field.k] || 0} onChange={(e) => handleFluidChange(f.id, field.k as keyof Fluid, parseFloat(e.target.value))} className="w-full p-2 border border-slate-300 rounded text-sm"/>
                         <span className="absolute right-2 top-2 text-xs text-slate-400">{field.u}</span>
                       </div>
                     </div>
                   ))}
                </div>
              </div>
            ))}
          </div>
        </div>
        <div className="p-5 border-t border-slate-200 bg-white flex justify-between">
          <button onClick={onReset} className="flex items-center gap-2 text-sm font-medium text-slate-600 hover:text-red-600"><RotateCcw size={16} /> Reset Defaults</button>
          <div className="flex gap-3">
            <button onClick={onClose} className="px-5 py-2 text-sm font-medium text-slate-600 hover:bg-slate-100 rounded-lg">Cancel</button>
            <button onClick={() => onSave(localFluids)} className="flex items-center gap-2 px-6 py-2 text-sm font-bold text-white bg-blue-600 hover:bg-blue-700 rounded-lg"><Save size={18} /> Save Database</button>
          </div>
        </div>
      </div>
    </div>
  );
};

/* ==============================================================================================
   SECTION 5: MAIN DASHBOARD
   ============================================================================================== */

const Dashboard = () => {
  // --- State ---
  const [simMode, setSimMode] = useState<'design' | 'rating'>('rating'); 
  const [inputs, setInputs] = useState({
    Tj: 85,             // Max Junction Temp Target (Design Mode)
    targetPower: 1000,  // Target Power Dissipation (Rating Mode)
    T_inlet: 40,        // Bulk Fluid Temp
    bypassFactor: 0.1,  // Leakage/Bypass
  });
  const [hsGeometry, setHsGeometry] = useState<Geometry>(INITIAL_GEOMETRY);
  const [fluids, setFluids] = useState<Fluid[]>(DEFAULT_FLUIDS);
  const [selectedFluids, setSelectedFluids] = useState<string[]>(['mpao', 'fuchs', 'dcf281']);
  
  // Modals
  const [isFluidModalOpen, setIsFluidModalOpen] = useState(false);
  const [isMethodModalOpen, setIsMethodModalOpen] = useState(false);
  
  // Data
  const [simData, setSimData] = useState<any[]>([]);
  const [rankingData, setRankingData] = useState<any[]>([]);
  const [viscosityData, setViscosityData] = useState<any[]>([]);
  const [densityData, setDensityData] = useState<any[]>([]);
  const [rayleighData, setRayleighData] = useState<any[]>([]);

  // --- Derived Metrics ---
  const derivedGeom = useMemo(() => {
    // 1. Fin Count (N)
    const N = Math.floor(hsGeometry.W / (hsGeometry.s + hsGeometry.t_fin));
    
    // 2. Surface Areas
    const A_fin = N * (2 * hsGeometry.h_fin * hsGeometry.L + hsGeometry.t_fin * hsGeometry.L);
    const A_base = N * (hsGeometry.s * hsGeometry.L); 
    const A_total = A_fin + A_base;
    
    // 3. Flow Area & Hydraulic Diameter
    const A_c_total = N * (hsGeometry.s * hsGeometry.h_fin);
    // Dh = 4 * Area / Perimeter = 2s (approx for parallel plates)
    const Dh = (2 * hsGeometry.s * hsGeometry.h_fin) / (hsGeometry.s + hsGeometry.h_fin);

    // 4. Spreading Resistance (Simplified Yovanovich)
    const A_die = hsGeometry.die_W * hsGeometry.die_L;
    const A_hs_base = hsGeometry.W * hsGeometry.L;
    const epsilon = Math.sqrt(A_die / A_hs_base); 
    const R_spreading = (1 - epsilon) / (hsGeometry.k_base * Math.sqrt(A_die) * Math.sqrt(Math.PI));
    
    const R_internal = hsGeometry.R_TIM + hsGeometry.R_die + R_spreading;

    return { N, A_total, A_fin, A_base, A_c_total, R_internal, Dh, R_spreading, A_die };
  }, [hsGeometry]);

  // --- Simulation Loop ---
  useEffect(() => {
    // Generate height sweep data (0.5m to 3.0m)
    const heightPoints = [];
    for (let h = 0.5; h <= 3.0; h += 0.1) {
      const point: any = { height: Number(h.toFixed(1)) };
      
      selectedFluids.forEach(fluidId => {
        const fluid = fluids.find(f => f.id === fluidId);
        if (fluid) {
          const res = simMode === 'rating' 
            ? solveForTj(inputs.targetPower, fluid, h, inputs.T_inlet, hsGeometry, derivedGeom, inputs.bypassFactor)
            : solveOperatingPoint(fluid, h, inputs.Tj, inputs.T_inlet, hsGeometry, derivedGeom, inputs.bypassFactor);

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

    // Generate Ranking (at fixed height 1.5m)
    const rankH = 1.5;
    const ranks = selectedFluids.map(fluidId => {
      const fluid = fluids.find(f => f.id === fluidId);
      if (!fluid) return null;
      
      const res = simMode === 'rating'
        ? solveForTj(inputs.targetPower, fluid, rankH, inputs.T_inlet, hsGeometry, derivedGeom, inputs.bypassFactor)
        : solveOperatingPoint(fluid, rankH, inputs.Tj, inputs.T_inlet, hsGeometry, derivedGeom, inputs.bypassFactor);
        
      return res ? { ...res, name: fluid.name, color: fluid.color } : null;
    }).filter(Boolean).sort((a: any, b: any) => simMode === 'rating' ? a.T_j - b.T_j : b.Q - a.Q);
    
    setRankingData(ranks);

    // Generate Viscosity Data
    const viscPoints = [];
    for (let T = 30; T <= 120; T += 5) {
      const point: any = { temp: T };
      selectedFluids.forEach(id => {
        const fluid = fluids.find(f => f.id === id);
        if (fluid) point[id] = getViscosity(fluid, T) * 1e6; // to cSt
      });
      viscPoints.push(point);
    }
    setViscosityData(viscPoints);

    // Generate Density Data
    const densPoints = [];
    for (let T = 30; T <= 120; T += 5) {
      const point: any = { temp: T };
      selectedFluids.forEach(id => {
        const fluid = fluids.find(f => f.id === id);
        if (fluid) {
          // Calculate density using linear expansion from 40C reference (consistent with solver)
          point[id] = fluid.rho * (1 - fluid.beta * (T - 40));
        }
      });
      densPoints.push(point);
    }
    setDensityData(densPoints);

    // Generate Rayleigh Data (Ra)
    const raPoints = [];
    for (let T = 30; T <= 120; T += 5) {
      const point: any = { temp: T };
      selectedFluids.forEach(id => {
        const fluid = fluids.find(f => f.id === id);
        if (fluid) {
          const nu = getViscosity(fluid, T);
          const rho = fluid.rho * (1 - fluid.beta * (T - 40));
          const cp = fluid.cp * (1 + 0.0015 * (T - 40));
          const k = fluid.k; 
          
          const alpha = k / (rho * cp);
          
          // Ra_s based on fin spacing (s) and fixed dT=50C for comparison
          const dT = 50; 
          const s = hsGeometry.s;
          // Ra = (g * beta * dT * s^3) / (nu * alpha)
          const Ra = (CONSTANTS.G * fluid.beta * dT * Math.pow(s, 3)) / (nu * alpha);
          
          point[id] = Ra;
        }
      });
      raPoints.push(point);
    }
    setRayleighData(raPoints);

  }, [inputs, selectedFluids, hsGeometry, derivedGeom, fluids, simMode]);

  // --- Helpers ---
  const toggleFluid = (id: string) => setSelectedFluids(prev => prev.includes(id) ? prev.filter(f => f !== id) : [...prev, id]);

  const downloadCSV = () => {
    let csv = "Height(m),";
    selectedFluids.forEach(f => csv += `${f} Power(W),${f} Tj(C),${f} Tout(C),${f} Velocity(m/s),`);
    csv += "\n";
    simData.forEach(row => {
      csv += `${row.height},`;
      selectedFluids.forEach(f => csv += `${row[`${f}_Q`]?.toFixed(2)},${row[`${f}_Tj`]?.toFixed(2)},${row[`${f}_Tout`]?.toFixed(2)},${row[`${f}_u`]?.toFixed(4)},`);
      csv += "\n";
    });
    const a = document.createElement('a');
    a.href = window.URL.createObjectURL(new Blob([csv], { type: 'text/csv' }));
    a.download = `simulation_${simMode}.csv`;
    a.click();
  };

  return (
    <div className="min-h-screen bg-slate-50 font-sans text-slate-900">
      
      <FluidEditModal 
        isOpen={isFluidModalOpen} onClose={() => setIsFluidModalOpen(false)}
        fluids={fluids} onSave={(newFluids: Fluid[]) => { setFluids(newFluids); setIsFluidModalOpen(false); }}
        onReset={() => { setFluids(DEFAULT_FLUIDS); setIsFluidModalOpen(false); }}
      />
      
      {/* Methodology Modal omitted for brevity, logic remains identical */}

      {/* --- Header --- */}
      <header className="bg-slate-900 text-white p-6 shadow-lg">
        <div className="max-w-7xl mx-auto flex justify-between items-center">
          <div>
            <h1 className="text-2xl font-bold flex items-center gap-3"><Wind className="text-blue-400" /> Natural Convection Simulator</h1>
            <p className="text-slate-400 text-sm mt-1">Single-Phase Immersion Cooling • {simMode === 'rating' ? 'Fixed Power Mode' : 'Fixed Temp Mode'}</p>
          </div>
          <div className="flex items-center gap-6">
             <div className="text-right hidden md:block pl-6 border-l border-slate-700">
               {simMode === 'rating' ? (
                 <> <div className="text-xs text-slate-400">Target Power</div> <div className="text-xl font-bold text-green-400">{inputs.targetPower} W</div> </>
               ) : (
                 <> <div className="text-xs text-slate-400">Max Temp</div> <div className="text-xl font-bold text-red-400">{inputs.Tj} °C</div> </>
               )}
            </div>
          </div>
        </div>
      </header>

      {/* --- Main Grid --- */}
      <main className="max-w-7xl mx-auto p-6 grid grid-cols-1 lg:grid-cols-12 gap-6">
        
        {/* Left Column: Controls */}
        <div className="lg:col-span-3 space-y-6">
          
          {/* 1. Operating Point */}
          <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5">
            <h3 className="font-semibold text-slate-800 flex items-center gap-2 mb-4"><Settings size={18} /> Operating Conditions</h3>
            <div className="flex bg-slate-100 p-1 rounded-lg mb-6">
              <button onClick={() => setSimMode('rating')} className={`flex-1 py-1.5 text-xs font-semibold rounded-md transition-all ${simMode === 'rating' ? 'bg-white text-blue-600 shadow-sm' : 'text-slate-500 hover:text-slate-700'}`}>Fixed Power</button>
              <button onClick={() => setSimMode('design')} className={`flex-1 py-1.5 text-xs font-semibold rounded-md transition-all ${simMode === 'design' ? 'bg-white text-blue-600 shadow-sm' : 'text-slate-500 hover:text-slate-700'}`}>Fixed Temp</button>
            </div>
            <div className="space-y-6">
              {simMode === 'design' ? (
                <div>
                  <label className="flex justify-between text-sm font-medium text-slate-700 mb-2">Max Junction Temp <span className="text-blue-600">{inputs.Tj}°C</span></label>
                  <input type="range" min="60" max="120" step="1" value={inputs.Tj} onChange={(e) => setInputs({...inputs, Tj: Number(e.target.value)})} className="w-full h-2 bg-slate-200 rounded-lg cursor-pointer accent-blue-600"/>
                </div>
              ) : (
                <div>
                  <label className="flex justify-between text-sm font-medium text-slate-700 mb-2">Target Power <span className="text-blue-600">{inputs.targetPower} W</span></label>
                  <input type="range" min="500" max="2500" step="50" value={inputs.targetPower} onChange={(e) => setInputs({...inputs, targetPower: Number(e.target.value)})} className="w-full h-2 bg-slate-200 rounded-lg cursor-pointer accent-blue-600"/>
                </div>
              )}
              <div>
                <label className="flex justify-between text-sm font-medium text-slate-700 mb-2">Inlet Fluid Temp <span className="text-blue-600">{inputs.T_inlet}°C</span></label>
                <input type="range" min="20" max="60" step="1" value={inputs.T_inlet} onChange={(e) => setInputs({...inputs, T_inlet: Number(e.target.value)})} className="w-full h-2 bg-slate-200 rounded-lg cursor-pointer accent-blue-600"/>
              </div>
              <div>
                <label className="flex justify-between text-sm font-medium text-slate-700 mb-2">Ducting Bypass <span className="text-orange-600">{(inputs.bypassFactor * 100).toFixed(0)}%</span></label>
                <input type="range" min="0" max="0.5" step="0.05" value={inputs.bypassFactor} onChange={(e) => setInputs({...inputs, bypassFactor: Number(e.target.value)})} className="w-full h-2 bg-slate-200 rounded-lg cursor-pointer accent-orange-500"/>
              </div>
            </div>
          </div>

          {/* 2. Geometry */}
          <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5">
            <h3 className="font-semibold text-slate-800 flex items-center gap-2 mb-4"><Ruler size={18} /> Heat Sink Geometry</h3>
            <div className="space-y-1">
               <InputField label="Height (Vertical)" value={hsGeometry.L} unit="m" step={0.01} onChange={(v: number) => setHsGeometry(p => ({...p, L: v}))}/>
               <InputField label="Width (W)" value={hsGeometry.W} unit="m" step={0.005} onChange={(v: number) => setHsGeometry(p => ({...p, W: v}))}/>
               <InputField label="Fin Height" value={hsGeometry.h_fin} unit="m" step={0.001} onChange={(v: number) => setHsGeometry(p => ({...p, h_fin: v}))}/>
               <InputField label="Fin Spacing" value={hsGeometry.s} unit="m" step={0.0005} onChange={(v: number) => setHsGeometry(p => ({...p, s: v}))}/>
               <InputField label="Fin Thickness" value={hsGeometry.t_fin} unit="m" step={0.0001} onChange={(v: number) => setHsGeometry(p => ({...p, t_fin: v}))}/>
            </div>
            <div className="mt-4 pt-4 border-t border-slate-100 space-y-2">
               <h4 className="text-xs font-bold text-slate-500 uppercase flex items-center gap-1"><Cpu size={12}/> Die Geometry</h4>
               <div className="grid grid-cols-2 gap-2">
                  <InputField label="Die Width" value={hsGeometry.die_W} unit="m" step={0.001} onChange={(v: number) => setHsGeometry(p => ({...p, die_W: v}))}/>
                  <InputField label="Die Height" value={hsGeometry.die_L} unit="m" step={0.001} onChange={(v: number) => setHsGeometry(p => ({...p, die_L: v}))}/>
               </div>
            </div>
            <div className="mt-4 pt-4 border-t border-slate-100 space-y-2">
               <h4 className="text-xs font-bold text-slate-500 uppercase">Derived Metrics</h4>
               <div className="flex justify-between text-sm"><span className="text-slate-600">Hydraulic Dia:</span><span className="font-mono font-bold">{(derivedGeom.Dh * 1000).toFixed(2)} mm</span></div>
               <div className="flex justify-between text-sm"><span className="text-slate-600">Spreading Res:</span><span className="font-mono font-bold text-orange-600">{(derivedGeom.R_spreading).toFixed(4)} K/W</span></div>
            </div>
          </div>

          {/* 3. Fluid Select */}
          <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5">
            <div className="flex justify-between items-center mb-4">
              <h3 className="font-semibold text-slate-800 flex items-center gap-2"><Droplets size={18} /> Fluids</h3>
              <button onClick={() => setIsFluidModalOpen(true)} className="text-xs bg-blue-50 hover:bg-blue-100 text-blue-700 px-2 py-1 rounded border border-blue-200 transition-colors flex items-center gap-1"><Edit2 size={12} /> Edit</button>
            </div>
            <div className="space-y-2">
              {fluids.map(fluid => (
                <div key={fluid.id} onClick={() => toggleFluid(fluid.id)} className={`flex items-center justify-between p-2 rounded-md cursor-pointer transition-colors border ${selectedFluids.includes(fluid.id) ? 'bg-blue-50 border-blue-200' : 'hover:bg-slate-50 border-transparent'}`}>
                  <div className="flex items-center gap-2"><div className="w-3 h-3 rounded-full" style={{backgroundColor: fluid.color}}></div><span className="text-sm font-medium text-slate-700">{fluid.name}</span></div>
                  {selectedFluids.includes(fluid.id) && <div className="w-2 h-2 bg-blue-500 rounded-full"></div>}
                </div>
              ))}
            </div>
          </div>
        </div>

        {/* Right Column: Charts */}
        <div className="lg:col-span-9 space-y-6">
          
          {/* Main Chart */}
          <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
            <div className="flex justify-between items-center mb-6">
              <div>
                <h2 className="text-lg font-bold text-slate-900">{simMode === 'rating' ? 'Junction Temperature vs Tank Height' : 'Power Dissipation vs Tank Height'}</h2>
                <p className="text-sm text-slate-500">{simMode === 'rating' ? `Solving for Tj at ${inputs.targetPower}W load.` : `Solving for Max Power at ${inputs.Tj}°C limit.`}</p>
              </div>
              <button onClick={downloadCSV} className="flex items-center gap-2 px-3 py-1 bg-slate-100 hover:bg-slate-200 text-slate-600 text-xs font-bold rounded-lg transition-colors"><Download size={14}/> CSV</button>
            </div>
            <div className="h-[350px] w-full">
              <ResponsiveContainer width="100%" height="100%">
                <ComposedChart data={simData} margin={{ top: 5, right: 30, left: 20, bottom: 25 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                  <XAxis dataKey="height" label={{ value: 'Tank Height (m)', position: 'insideBottom', offset: -15 }} type="number" domain={[0.5, 3.0]} tickCount={6}/>
                  <YAxis label={{ value: simMode === 'rating' ? 'Temp (°C)' : 'Power (W)', angle: -90, position: 'insideLeft' }} domain={['auto', 'auto']}/>
                  <Tooltip contentStyle={{ backgroundColor: '#fff', borderRadius: '8px', border: '1px solid #e2e8f0' }} labelFormatter={(label) => `Height: ${label} m`} formatter={(value: number, name: string) => [simMode === 'rating' ? `${value.toFixed(1)} °C` : `${value.toFixed(0)} W`, name.split('_')[0]]}/>
                  <Legend verticalAlign="top" height={36}/>
                  {simMode === 'rating' ? <ReferenceLine y={108} stroke="red" strokeDasharray="3 3" label={{ position: 'insideRight', value: 'Max 108°C' }} /> : <ReferenceLine y={inputs.targetPower} stroke="green" strokeDasharray="3 3" label="Target" />}
                  {fluids.map(fluid => selectedFluids.includes(fluid.id) && <Line key={fluid.id} type="monotone" dataKey={simMode === 'rating' ? `${fluid.id}_Tj` : `${fluid.id}_Q`} name={fluid.name} stroke={fluid.color} strokeWidth={3} dot={false}/>)}
                </ComposedChart>
              </ResponsiveContainer>
            </div>
          </div>

          {/* Secondary Charts */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
              <h2 className="text-lg font-bold text-slate-900 mb-1">Induced Velocity</h2>
              <p className="text-sm text-slate-500 mb-4">Net velocity through fins (after leakage)</p>
              <div className="h-[250px] w-full">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={simData}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="height" />
                    <YAxis label={{ value: 'm/s', angle: -90, position: 'insideLeft' }} />
                    <Tooltip formatter={(val: number) => val.toFixed(3)} />
                    {fluids.map(fluid => selectedFluids.includes(fluid.id) && <Line key={fluid.id} type="monotone" dataKey={`${fluid.id}_u`} stroke={fluid.color} dot={false} strokeWidth={2}/>)}
                  </LineChart>
                </ResponsiveContainer>
              </div>
            </div>
            <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
              <h2 className="text-lg font-bold text-slate-900 mb-1">Viscosity Profile</h2>
              <p className="text-sm text-slate-500 mb-4">ASTM D341 Model (Log Scale)</p>
              <div className="h-[250px] w-full">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={viscosityData}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="temp" label={{ value: '°C', position: 'insideBottom', offset: -5 }}/>
                    <YAxis scale="log" domain={['auto', 'auto']} allowDataOverflow label={{ value: 'cSt', angle: -90, position: 'insideLeft' }} />
                    <Tooltip formatter={(val: number) => val.toFixed(2)} />
                    {fluids.map(fluid => selectedFluids.includes(fluid.id) && <Line key={fluid.id} type="monotone" dataKey={fluid.id} stroke={fluid.color} dot={false} strokeWidth={2}/>)}
                  </LineChart>
                </ResponsiveContainer>
              </div>
            </div>
            <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
              <h2 className="text-lg font-bold text-slate-900 mb-1">Density Profile</h2>
              <p className="text-sm text-slate-500 mb-4">Linear Expansion (Ref 40°C)</p>
              <div className="h-[250px] w-full">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={densityData}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="temp" label={{ value: '°C', position: 'insideBottom', offset: -5 }}/>
                    <YAxis domain={['auto', 'auto']} label={{ value: 'kg/m³', angle: -90, position: 'insideLeft' }} />
                    <Tooltip formatter={(val: number) => val.toFixed(1)} />
                    {fluids.map(fluid => selectedFluids.includes(fluid.id) && <Line key={fluid.id} type="monotone" dataKey={fluid.id} stroke={fluid.color} dot={false} strokeWidth={2}/>)}
                  </LineChart>
                </ResponsiveContainer>
              </div>
            </div>
             <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
              <h2 className="text-lg font-bold text-slate-900 mb-1">Rayleigh Number</h2>
              <p className="text-sm text-slate-500 mb-4">Dimensionless (Based on Spacing, ΔT=50°C)</p>
              <div className="h-[250px] w-full">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={rayleighData}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="temp" label={{ value: '°C', position: 'insideBottom', offset: -5 }}/>
                    <YAxis domain={['auto', 'auto']} tickFormatter={(v) => v.toExponential(0)} label={{ value: 'Ra', angle: -90, position: 'insideLeft' }} />
                    <Tooltip formatter={(val: number) => val.toExponential(2)} />
                    {fluids.map(fluid => selectedFluids.includes(fluid.id) && <Line key={fluid.id} type="monotone" dataKey={fluid.id} stroke={fluid.color} dot={false} strokeWidth={2}/>)}
                  </LineChart>
                </ResponsiveContainer>
              </div>
            </div>
          </div>
          
          {/* Ranking Table */}
          <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
            <div className="p-6 border-b border-slate-200">
              <h2 className="text-lg font-bold text-slate-900">Performance Snapshot (Height = 1.5m)</h2>
            </div>
            <div className="overflow-x-auto">
              <table className="w-full text-left text-sm">
                <thead className="bg-slate-50 text-slate-600 font-semibold">
                  <tr>
                    <th className="p-4">Rank</th>
                    <th className="p-4">Fluid</th>
                    <th className="p-4">T_junc (°C)</th>
                    <th className="p-4">Power (W)</th>
                    <th className="p-4">T_out (°C)</th>
                    <th className="p-4">Vel (m/s)</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-100">
                  {rankingData.map((row, idx) => (
                    <tr key={idx} className="hover:bg-slate-50">
                      <td className="p-4 font-mono text-slate-400">#{idx + 1}</td>
                      <td className="p-4 font-medium flex items-center gap-2"><div className="w-3 h-3 rounded-full" style={{backgroundColor: row.color}}></div>{row.name}</td>
                      <td className="p-4 font-bold">{row.T_j.toFixed(1)}</td>
                      <td className="p-4 font-bold">{row.Q.toFixed(0)}</td>
                      <td className="p-4">{row.T_out.toFixed(1)}</td>
                      <td className="p-4">{row.u_ind.toFixed(3)}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>

        </div>
      </main>
    </div>
  );
};

const root = createRoot(document.getElementById('root')!);
root.render(<Dashboard />);