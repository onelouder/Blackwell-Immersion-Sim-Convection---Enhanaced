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
import { Settings, Thermometer, Wind, Activity, Info, Droplets, ArrowRight, Ruler, Calculator, Sliders, Edit2, X, Save, RotateCcw } from 'lucide-react';

// --- 1. Physics Constants & Fluid Database ---

const G = 9.81; // Gravity (m/s^2)

// Default Geometry defined in "Final Parameter Summary" PDF Page 8
const INITIAL_GEOMETRY = {
  L: 0.15, // Height of HS in vertical orientation (m)
  W: 0.08, // Width (m)
  s: 0.002, // Fin spacing (m)
  t_fin: 0.0005, // Fin thickness (m)
  h_fin: 0.02, // Fin height/length (m)
  base_thick: 0.004, // Base thickness (estimated)
  
  // Thermal Resistances (K/W)
  R_JC: 0.000390625, // Junction to Case
  R_CS: 0.000390625, // Case to Sink
  R_cond: 0.004667, // Conduction resistance
  
  // Fin Material (Copper)
  k_fin: 400, // W/m-K
};

// Fluid Properties
// nu_40 and nu_100 are Kinematic Viscosity in cSt (mm2/s) used for ASTM interpolation
const DEFAULT_FLUIDS = [
  { 
    id: 'dcf281', 
    name: 'DCF-281 (Novvi)', 
    rho: 793.36, 
    nu: 10.41e-6, // m2/s @ 40C
    nu_cSt_40: 10.41,
    nu_cSt_100: 2.81, // Updated by user request
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
  // If we don't have 100C data (custom fluids), fallback to constant
  if (!fluid.nu_cSt_100 || !fluid.nu_cSt_40) return fluid.nu;

  const T_K = T_C + 273.15;
  const T1 = 40 + 273.15;
  const T2 = 100 + 273.15;

  const v1 = fluid.nu_cSt_40;
  const v2 = fluid.nu_cSt_100;

  // ASTM D341: log(log(v + 0.7)) = A - B * log(T)
  // Z = log10(log10(v + 0.7))
  const Z1 = Math.log10(Math.log10(v1 + 0.7));
  const Z2 = Math.log10(Math.log10(v2 + 0.7));

  const B = (Z1 - Z2) / (Math.log10(T2) - Math.log10(T1));
  const A = Z1 + B * Math.log10(T1);

  const Z_target = A - B * Math.log10(T_K);
  // v = 10^(10^Z) - 0.7
  const v_cSt = Math.pow(10, Math.pow(10, Z_target)) - 0.7;

  // Return in m^2/s
  return Math.max(0.1, v_cSt) * 1e-6; // Clamp min viscosity to avoid singularities
};


// Developing Laminar Flow Correlation (Hausen / Sieder-Tate approx)
// Includes Sieder-Tate Correction for property variation
const calculateLaminarNu = (Re: number, Pr: number, Dh: number, L: number, mu_bulk: number, mu_wall: number) => {
  // Nu = Nu_fd + (0.0668 * (Dh/L) * Re * Pr) / (1 + 0.04 * [(Dh/L) * Re * Pr]^(2/3))
  // Sieder-Tate Factor: (mu_b / mu_w)^0.14
  
  const Nu_fd = 7.541; // Rectangular duct aspect ratio approx
  const Gz = (Dh / L) * Re * Pr; // Graetz number approx term
  
  const numerator = 0.0668 * Gz;
  const denominator = 1 + 0.04 * Math.pow(Gz, 2/3);
  
  const Nu_iso = Nu_fd + (numerator / denominator);
  
  // Correction for hot wall / cool fluid (Heating mode)
  // Wall is hotter, so mu_wall < mu_bulk. Ratio > 1. Nu increases.
  const visc_correction = Math.pow(mu_bulk / mu_wall, 0.14);
  
  return Nu_iso * visc_correction;
};

// Fin Efficiency
const calculateFinEff = (h: number, geom: typeof INITIAL_GEOMETRY) => {
  // m = sqrt(2h / (k_fin * t_fin))
  const m = Math.sqrt((2 * h) / (geom.k_fin * geom.t_fin));
  const mh = m * geom.h_fin;
  // tanh(mh) / mh
  if (mh === 0) return 1;
  return Math.tanh(mh) / mh;
};

// Iterative Solver for a single Operating Point
const solveOperatingPoint = (
  fluid: any, 
  H_tank: number, 
  Tj: number, 
  T_inlet: number,
  geom: typeof INITIAL_GEOMETRY,
  derived: { N: number, A_total: number, A_fin: number, A_base: number, A_c_total: number, R_internal: number, Dh: number }
) => {
  // Constant properties
  const beta = fluid.beta; 

  // Iteration State
  let T_wall = Tj - 2; // Initial guess
  let T_out = T_inlet + 5; // Initial guess
  
  // Results
  let Q_dissipated = 0;
  let u_induced = 0;
  let h_conv = 0;
  let Re = 0;
  let Ra = 0;
  let nu_op = fluid.nu; // Operating viscosity for reporting
  
  let converged = false;
  let iterations = 0;
  
  // Relaxation factor
  const alpha = 0.2; 
  
  // Pressure Loss Coefficients
  // Shah & London Developing Flow Friction Factor
  // f_app * Re = 3.44 / sqrt(L+) ... simplified approx
  const calcFrictionFactor = (Re_local: number, Dh: number, L: number) => {
    if (Re_local < 1) return 24;
    // Dimensionless length L+ = L / (Dh * Re)
    const L_plus = L / (Dh * Re_local);
    // Shah & London approx for parallel plates
    // f_Re = 24 + (0.674 / (4 * L_plus)) / (1 + 0.00029 / (4 * L_plus)^2) ... simplified curve fit
    // Using a simpler composite fit:
    // Fully developed f*Re = 24
    // Developing region adds penalty
    const fRe_fd = 24;
    const fRe_dev = 3.44 / Math.sqrt(L_plus);
    const fRe_app = Math.sqrt(fRe_fd*fRe_fd + fRe_dev*fRe_dev);
    
    return fRe_app / Re_local;
  };
  
  // Contraction and Expansion Losses (Kc, Ke) based on Area Ratio (sigma)
  // sigma = flow_area / frontal_area
  // For heat sink in tank, frontal area is W * (h_fin + base). This is rough approx.
  // Better: sigma = A_c_total / (W * (h_fin + t_fin))
  // Assuming 'W' is total width and fins fill it.
  const sigma = (geom.s * derived.N) / (geom.W); 
  // Kays & London approx for laminar flow
  const Kc = 0.42 * (1 - sigma * sigma);
  const Ke = (1 - sigma) * (1 - sigma);
  const K_total = Kc + Ke;

  while (!converged && iterations < 50) {
    // A. Property Evaluation
    const T_mean_fluid = (T_inlet + T_out) / 2;
    
    // 1. Calculate Viscosity at Operating Temps
    const nu_bulk = getViscosity(fluid, T_mean_fluid);
    const nu_wall = getViscosity(fluid, T_wall);
    nu_op = nu_bulk; // Store for reporting

    // 2. Calculate Density Correction (Linear Boussinesq)
    // rho_bulk = rho_40 * (1 - beta * (T_mean - 40))
    const rho_bulk = fluid.rho * (1 - beta * (T_mean_fluid - 40));
    
    // B. Induced Velocity (Quadratic Solver)
    // IMPROVEMENT: Split Buoyancy Integral
    // Core (Heat Sink): T_mean - T_inlet over height L
    // Chimney (Plume): T_out - T_inlet over height (H_tank - L)
    
    const dT_core = Math.max(0.1, T_mean_fluid - T_inlet);
    const dT_chimney = Math.max(0.1, T_out - T_inlet);
    
    const H_chimney = Math.max(0, H_tank - geom.L);
    
    const dP_buoy_core = rho_bulk * G * beta * dT_core * geom.L;
    const dP_buoy_chimney = rho_bulk * G * beta * dT_chimney * H_chimney;
    const dP_buoy_total = dP_buoy_core + dP_buoy_chimney;

    // Resisting Pressure:
    // 1. Viscous: f_corrected * (1/2 rho u^2) * (L/Dh) * 4 ? No, standard D-W is f * (L/D) * dyn_head
    // or Hagen-Poiseuille equiv. 
    // Using Darcy-Weisbach: dP = f * (L/Dh) * (1/2 rho u^2)
    // where f is Apparent Friction Factor for developing flow
    
    // Iterative velocity required for Re-dependent Friction
    // Simplified: Assume laminar relationship first to get coefficients
    // But f depends on u (via Re). 
    // Let's use the quadratic form derived from f = C/Re
    // dP_visc = (C/Re) * (L/Dh) * 0.5 * rho * u^2 = (C * nu / (u * Dh)) * (L/Dh) * 0.5 * rho * u^2
    // = C * nu * L * 0.5 * rho * u / Dh^2
    // Linearly proportional to u for fully developed. 
    // However, developing flow has sqrt terms. 
    // Let's iterate u specifically or use a linearized friction coefficient at current Re guess.
    
    const u_guess = u_induced > 0 ? u_induced : 0.1;
    const Re_guess = (u_guess * derived.Dh) / nu_bulk;
    const f_app = calcFrictionFactor(Re_guess, derived.Dh, geom.L);
    
    // Viscosity Correction for Friction (Liquid Heating)
    // Wall is hotter -> lower visc near wall -> lower friction than iso
    const mu_bulk = nu_bulk * rho_bulk;
    const mu_wall_calc = nu_wall * rho_bulk;
    const visc_ratio = mu_wall_calc / mu_bulk;
    const friction_visc_correction = Math.pow(visc_ratio, 0.58); 
    
    const f_effective = f_app * friction_visc_correction;

    // dP_loss = (K + f*L/Dh) * 0.5 * rho * u^2
    const K_flow = K_total + (f_effective * geom.L / derived.Dh);
    
    // Balance: dP_buoy = K_flow * 0.5 * rho * u^2
    // u = sqrt( 2 * dP_buoy / (rho * K_flow) )
    const u_new = Math.sqrt( (2 * dP_buoy_total) / (rho_bulk * K_flow) );
    
    u_induced = u_induced * 0.5 + u_new * 0.5; // Smooth update
    if (u_induced < 0.0001) u_induced = 0.0001;
    
    // C. Heat Transfer
    Re = (u_induced * derived.Dh) / nu_bulk;
    
    // Calculate Rayleigh Number for reporting
    // Ra = Gr * Pr
    // Gr = g * beta * (Tw - Tin) * L^3 / nu^2
    const Gr = (G * beta * (T_wall - T_inlet) * Math.pow(geom.L, 3)) / (nu_bulk * nu_bulk);
    const Pr = (mu_bulk * fluid.cp) / fluid.k;
    Ra = Gr * Pr;
    
    // Calculate Nu with Sieder-Tate Correction (mu_bulk / mu_wall)
    // Note: Passing dynamic mu values
    const Nu = calculateLaminarNu(Re, Pr, derived.Dh, geom.L, mu_bulk, mu_wall_calc);
    
    const h_bare = (Nu * fluid.k) / derived.Dh;
    
    // Fin Efficiency
    const eta_fin = calculateFinEff(h_bare, geom);
    const eta_o = 1 - (derived.A_fin / derived.A_total) * (1 - eta_fin);
    const h_eff = h_bare * eta_o;
    h_conv = h_eff;

    // D. Energy Balance
    // Q_conv = h * A * (T_wall - T_mean_fluid)
    const dT_conv = Math.max(0.1, T_wall - T_mean_fluid);
    const Q_conv = h_eff * derived.A_total * dT_conv;

    // E. Outlet Temp & Wall Temp Calculation
    const m_dot = rho_bulk * u_induced * derived.A_c_total;
    
    let T_out_calc = T_inlet;
    if (m_dot > 0) {
      T_out_calc = T_inlet + (Q_conv / (m_dot * fluid.cp));
    }

    const T_wall_calc = Tj - Q_conv * derived.R_internal;

    // F. Convergence
    const err_Tout = Math.abs(T_out_calc - T_out);
    const err_Twall = Math.abs(T_wall_calc - T_wall);

    if (err_Tout < 0.01 && err_Twall < 0.01) {
      converged = true;
      Q_dissipated = Q_conv;
      T_out = T_out_calc;
      T_wall = T_wall_calc;
    } else {
      T_out = T_out * (1 - alpha) + T_out_calc * alpha;
      T_wall = T_wall * (1 - alpha) + T_wall_calc * alpha;
      
      // Clamps
      if (T_out > T_wall - 0.5) T_out = T_wall - 0.5;
      if (T_wall > Tj) T_wall = Tj;
      if (T_wall < T_inlet) T_wall = T_inlet + 1;
    }
    iterations++;
  }

  return {
    Q: Q_dissipated,
    u_ind: u_induced,
    T_out: T_out,
    T_wall: T_wall,
    h_conv: h_conv,
    Re: Re,
    Ra: Ra, // Return Rayleigh
    nu_op: nu_op
  };
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
        className="w-full px-3 py-2 bg-white border border-slate-300 rounded text-sm text-slate-800 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-transparent"
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

  const handleFluidChange = (id: string, field: string, val: number) => {
    setLocalFluids((prev: any[]) => prev.map(f => {
      if (f.id === id) {
        return { ...f, [field]: val };
      }
      return f;
    }));
  };

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/50 backdrop-blur-sm p-4 animate-in fade-in duration-200">
      <div className="bg-white rounded-xl shadow-2xl w-full max-w-5xl max-h-[90vh] flex flex-col overflow-hidden">
        <div className="p-5 border-b border-slate-200 bg-slate-50 flex justify-between items-center">
          <div>
            <h2 className="text-xl font-bold text-slate-800 flex items-center gap-2">
              <Edit2 size={20} className="text-blue-600"/> Edit Fluid Properties
            </h2>
            <p className="text-sm text-slate-500 mt-1">Customize physical parameters for simulation scenarios.</p>
          </div>
          <button onClick={onClose} className="p-2 hover:bg-slate-200 rounded-full transition-colors text-slate-500">
            <X size={24} />
          </button>
        </div>

        <div className="flex-1 overflow-y-auto p-6 bg-slate-50/50">
          <div className="grid grid-cols-1 gap-4">
            {localFluids.map((f: any) => (
              <div key={f.id} className="bg-white border border-slate-200 rounded-lg p-4 shadow-sm hover:shadow-md transition-shadow">
                <div className="flex items-center gap-3 mb-4 pb-3 border-b border-slate-100">
                  <div className="w-4 h-4 rounded-full shadow-sm" style={{backgroundColor: f.color}}></div>
                  <h3 className="font-bold text-slate-800 text-lg">{f.name}</h3>
                  <span className="text-xs text-slate-400 font-mono bg-slate-100 px-2 py-1 rounded">ID: {f.id}</span>
                </div>
                
                <div className="grid grid-cols-2 md:grid-cols-6 gap-4">
                   <div className="space-y-1">
                     <label className="text-xs font-semibold text-slate-500 uppercase">Density (ρ)</label>
                     <div className="relative">
                       <input 
                         type="number" step="0.1" 
                         value={f.rho}
                         onChange={(e) => handleFluidChange(f.id, 'rho', parseFloat(e.target.value))}
                         className="w-full p-2 border border-slate-300 rounded text-sm"
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
                         className="w-full p-2 border border-slate-300 rounded text-sm"
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
                         className="w-full p-2 border border-slate-300 rounded text-sm"
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
                         className="w-full p-2 border border-slate-300 rounded text-sm"
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
                         className="w-full p-2 border border-slate-300 rounded text-sm"
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
                         className="w-full p-2 border border-slate-300 rounded text-sm"
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
              <Save size={18} /> Save Changes
            </button>
          </div>
        </div>
      </div>
    </div>
  );
};

const Dashboard = () => {
  const [inputs, setInputs] = useState({
    Tj: 80, // deg C
    T_inlet: 40, // deg C
  });

  const [hsGeometry, setHsGeometry] = useState(INITIAL_GEOMETRY);
  const [fluids, setFluids] = useState(DEFAULT_FLUIDS);
  const [isFluidModalOpen, setIsFluidModalOpen] = useState(false);
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
    const R_internal = hsGeometry.R_JC + hsGeometry.R_CS + hsGeometry.R_cond;
    return { N, A_total, A_fin, A_base, A_c_total, R_internal, Dh };
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
    for (let h = 0.5; h <= 3.0; h += 0.1) {
      const point: any = { height: Number(h.toFixed(1)) };
      
      selectedFluids.forEach(fluidId => {
        const fluid = fluids.find(f => f.id === fluidId);
        if (fluid) {
          const res = solveOperatingPoint(fluid, h, inputs.Tj, inputs.T_inlet, hsGeometry, derivedGeom);
          point[`${fluidId}_Q`] = res.Q;
          point[`${fluidId}_u`] = res.u_ind;
          point[`${fluidId}_Tout`] = res.T_out;
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
      const res = solveOperatingPoint(fluid, rankH, inputs.Tj, inputs.T_inlet, hsGeometry, derivedGeom);
      return {
        name: fluid.name,
        color: fluid.color,
        Q: res.Q,
        T_out: res.T_out,
        u_ind: res.u_ind,
        Re: res.Re,
        Ra: res.Ra,
        nu_op: res.nu_op
      };
    }).filter(Boolean).sort((a: any, b: any) => b.Q - a.Q);
    setRankingData(ranks);

  }, [inputs, selectedFluids, hsGeometry, derivedGeom, fluids]);

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

  return (
    <div className="min-h-screen bg-slate-50 font-sans text-slate-900">
      
      <FluidEditModal 
        isOpen={isFluidModalOpen} 
        onClose={() => setIsFluidModalOpen(false)}
        fluids={fluids}
        onSave={handleFluidSave}
        onReset={handleFluidReset}
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
              Blackwell H200 Immersion Cooling • High Fidelity Physics • Variable Viscosity Model
            </p>
          </div>
          <div className="text-right hidden md:block">
            <div className="text-xs text-slate-400">Target Power</div>
            <div className="text-xl font-bold text-green-400">1200 W</div>
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
            
            <div className="space-y-6">
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
                <div className="flex justify-between text-xs text-slate-400 mt-1">
                  <span>60°C</span>
                  <span>120°C</span>
                </div>
              </div>

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
                <div className="flex justify-between text-xs text-slate-400 mt-1">
                  <span>20°C</span>
                  <span>50°C</span>
                </div>
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
            
            <div className="mt-4 pt-4 border-t border-slate-100 space-y-2">
               <h4 className="text-xs font-bold text-slate-500 uppercase">Derived Properties</h4>
               <div className="flex justify-between text-sm">
                 <span className="text-slate-600">Fin Count:</span>
                 <span className="font-mono font-bold">{derivedGeom.N}</span>
               </div>
               <div className="flex justify-between text-sm">
                 <span className="text-slate-600">Hydraulic Dia:</span>
                 <span className="font-mono font-bold">{(derivedGeom.Dh * 1000).toFixed(2)} mm</span>
               </div>
               <div className="flex justify-between text-sm">
                 <span className="text-slate-600">Total Area:</span>
                 <span className="font-mono font-bold">{derivedGeom.A_total.toFixed(3)} m²</span>
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
          
          {/* Chart 1: Power Dissipation */}
          <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
            <div className="flex justify-between items-center mb-6">
              <div>
                <h2 className="text-lg font-bold text-slate-900">Steady-State Power Dissipation (Q)</h2>
                <p className="text-sm text-slate-500">Calculated Q vs. Effective Tank Chimney Height (H_tank)</p>
              </div>
              <div className="flex gap-2">
                 <span className="px-3 py-1 bg-green-100 text-green-800 text-xs font-bold rounded-full">Goal: 1200W</span>
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
                    label={{ value: 'Power Q (W)', angle: -90, position: 'insideLeft' }} 
                    domain={[0, 'auto']}
                  />
                  <Tooltip 
                    contentStyle={{ backgroundColor: '#fff', borderRadius: '8px', border: '1px solid #e2e8f0', boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)' }}
                    labelFormatter={(label) => `Tank Height: ${label} m`}
                    formatter={(value: number, name: string) => [`${value.toFixed(0)} W`, name]}
                  />
                  <Legend verticalAlign="top" height={36}/>
                  
                  {/* Reference Lines */}
                  <ReferenceLine y={1200} stroke="red" strokeDasharray="3 3" label={{ position: 'insideRight', value: 'Target 1200W', fill: 'red', fontSize: 12 }} />
                  <ReferenceLine y={700} stroke="gray" strokeDasharray="3 3" label={{ position: 'insideRight', value: 'Min 700W', fill: 'gray', fontSize: 12 }} />

                  {fluids.map(fluid => (
                    selectedFluids.includes(fluid.id) && (
                      <Line 
                        key={fluid.id}
                        type="monotone" 
                        dataKey={`${fluid.id}_Q`} 
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
              <h2 className="text-lg font-bold text-slate-900 mb-1">Induced Velocity</h2>
              <p className="text-sm text-slate-500 mb-4">Chimney effect velocity vs Tank Height</p>
              <div className="h-[250px] w-full">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={simData} margin={{ top: 5, right: 10, left: 0, bottom: 20 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis dataKey="height" label={{ value: 'Tank Height (m)', position: 'insideBottom', offset: -10, fontSize: 12 }} />
                    <YAxis label={{ value: 'Velocity (m/s)', angle: -90, position: 'insideLeft', fontSize: 12 }} />
                    <Tooltip 
                      labelFormatter={(label) => `Tank Height: ${label} m`}
                      formatter={(val: number, name: string) => [`${val.toFixed(3)} m/s`, name]} 
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
                      formatter={(val: number, name: string) => [`${val.toFixed(1)} °C`, name]} 
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
              <p className="text-sm text-slate-500">Snapshot at Tank Height = 1.5m. Note: 'Op. Visc' shows thermal thinning effect.</p>
            </div>
            <div className="overflow-x-auto">
              <table className="w-full text-left text-sm">
                <thead className="bg-slate-50 text-slate-600 font-semibold">
                  <tr>
                    <th className="p-4">Rank</th>
                    <th className="p-4">Fluid</th>
                    <th className="p-4">Power (Q)</th>
                    <th className="p-4">Outlet T</th>
                    <th className="p-4">Velocity</th>
                    <th className="p-4">Re</th>
                    <th className="p-4">Rayleigh</th>
                    <th className="p-4">Op. Visc (cSt)</th>
                    <th className="p-4">Rec.</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-100">
                  {rankingData.map((row, idx) => (
                    <tr key={idx} className="hover:bg-slate-50">
                      <td className="p-4 font-mono text-slate-400">#{idx + 1}</td>
                      <td className="p-4 font-medium flex items-center gap-3">
                        <div className="w-3 h-3 rounded-full shadow-sm" style={{backgroundColor: row.color}}></div>
                        <span className="text-slate-900">{row.name}</span>
                      </td>
                      <td className="p-4 font-bold text-slate-800">{row.Q.toFixed(1)} W</td>
                      <td className="p-4">{row.T_out.toFixed(1)} °C</td>
                      <td className="p-4">{row.u_ind.toFixed(3)} m/s</td>
                      <td className="p-4 font-mono text-xs">{row.Re.toFixed(0)}</td>
                      <td className="p-4 font-mono text-xs text-slate-500">{(row.Ra).toExponential(2)}</td>
                      <td className="p-4 font-mono text-xs text-slate-500">{(row.nu_op * 1e6).toFixed(2)}</td>
                      <td className="p-4">
                        {row.Q >= 1200 ? (
                          <span className="inline-flex items-center px-2 py-1 rounded-full text-xs font-medium bg-green-100 text-green-800">
                            Excellent
                          </span>
                        ) : row.Q >= 1000 ? (
                          <span className="inline-flex items-center px-2 py-1 rounded-full text-xs font-medium bg-yellow-100 text-yellow-800">
                            Good
                          </span>
                        ) : (
                          <span className="inline-flex items-center px-2 py-1 rounded-full text-xs font-medium bg-red-100 text-red-800">
                            Insufficient
                          </span>
                        )}
                      </td>
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

// Mount React App
const root = createRoot(document.getElementById('root')!);
root.render(<Dashboard />);