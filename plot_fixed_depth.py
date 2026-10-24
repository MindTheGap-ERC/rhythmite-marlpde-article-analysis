"""Fixed-depth analysis of existing output; no model runs are performed.

VTU is ASCII space-time output from Rhythmite export_to_vtu2. Its depth
ordering is reversed and its horizontal coordinate is graphical, not time.
This reader loads one array at a time (several hundred MB peak memory),
checks dimensions/coordinates and compares each bottom series with ASCII.
Export rounding tolerance includes both ASCII (six significant digits) and
VTU (six decimal places). Run normally, not with Python -O: assertions are
intentional input-validation failures. Raw model data remain external.
"""
import argparse,mmap,re,json
from pathlib import Path
import numpy as np
import h5py
from utility import oscillation_onset,first_zero,time_first_minimum
from plot_rhythmite_times_against_length import plot_times

parser = argparse.ArgumentParser(description="Extract Figure 5 diagnostics at 500 cm, validating VTU bottom records against ASCII.")
parser.add_argument("--rhythmite-runs", type=Path, required=True,
                    help="Root containing resolution/200nnx and Larger_system")
parser.add_argument("--marlpde-results", type=Path, required=True)
parser.add_argument("--output-stem", type=Path, default=Path("figures/rhythmite_marlpde_onset_500cm"))
parser.add_argument("--rhythmite-tstar-years", type=float, default=13190.)
parser.add_argument("--rhythmite-xstar-cm", type=float, default=1319.)
args = parser.parse_args()
R, M = args.rhythmite_runs, args.marlpde_results
args.output_stem.parent.mkdir(parents=True, exist_ok=True)

def diag(t,c,u=None):
 out={}
 for k,fun,v in [('onset',lambda t,v:oscillation_onset(t,v,prominence=.01,min_extrema=6)[0],c),('zero',first_zero,u),('minimum',time_first_minimum,u)]:
  if v is None:continue
  try:out[k]=float(fun(t,v))
  except ValueError:out[k]=None
 return out
rows=[]
for L in [500,625,750,875,1000,1125,1250,1625]:
 d=R/'resolution/200nnx' if L==500 else R/'Larger_system'/('L_'+str(L)+('_longer_integration' if L in (1125,1250) else ''))
 nx=int(L*200/500); a=np.loadtxt(d/f'solution_x_{nx-1:06d}.ascii'); nt=len(a)
 with (d/'solution.vtu').open('rb') as f,mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ) as mm:
  npnt=int(re.search(rb'NumberOfPoints=\s*[\'\"]\s*(\d+)',mm[:500])[1]); assert npnt==nt*nx,(L,npnt,nt,nx)
  # Read coordinates only at depth-block starts; the export reverses field depths.
  start=mm.find(b'<DataArray');start=mm.find(b'>',start)+1;end=mm.find(b'</DataArray>',start)
  xyz=np.fromstring(mm[start:end].decode(),sep=' ').reshape(nx,nt,3)
  grid=(xyz[-1,0,1]-xyz[:,0,1])*args.rhythmite_xstar_cm
  assert abs(grid[0]-L)<.002
  # Use exact uniform grid after checking rounded exported coordinates.
  depths=np.linspace(L,0,nx);assert np.max(abs(grid-depths))<.002
  del xyz
  j=500/(L/(nx-1));lo=min(int(np.floor(j)),nx-2);hi=lo+1;w=j-lo
  arrays={};checks={}
  for field,col in [('c_co',4),('u',6)]:
   token=("Name='"+field+"'").encode();pos=mm.find(token);assert pos>=0
   start=mm.find(b'>',pos)+1;end=mm.find(b'</DataArray>',start)
   v=np.fromstring(mm[start:end].decode(),sep=' ').reshape(nx,nt)
   err=np.abs(v[0]-a[:,col]);tol=5.1e-6*np.abs(a[:,col])+5.1e-7
   checks[field]={'max_abs_error':float(err.max()),'within_export_rounding':bool(np.all(err<=tol))}
   assert checks[field]['within_export_rounding'],(L,field,checks[field])
   arrays[field]=(1-w)*v[nx-1-lo]+w*v[nx-1-hi]
   arrays[field+'_bottom']=v[0].copy()
   del v
 t=a[:,0]*args.rhythmite_tstar_years/1000
 with h5py.File(M/f'L_{L}/LMAHeureuxPorosityDiff.hdf5','r') as h:
  fs=json.loads(json.loads(h.attrs['field_attributes'])['fields']);g=json.loads(fs[0]['grid']);n=g['shape'][0];bounds=np.array(g['bounds'][0])*float(h.attrs['Xstar']);dx=(bounds[1]-bounds[0])/n
  j=(500-bounds[0])/dx-.5;ml=min(max(int(np.floor(j)),0),n-2);mw=j-ml
  mc=(1-mw)*h['data'][:,3,ml]+mw*h['data'][:,3,ml+1];mt=h['times'][:]*float(h.attrs['Tstar'])/1000
 result={'L':L,'rhythmite_vtu':str(d/'solution.vtu'),'rhythmite_ascii':str(d/f'solution_x_{nx-1:06d}.ascii'),'marlpde_hdf5':str(M/f'L_{L}/LMAHeureuxPorosityDiff.hdf5'),'rhythmite_nodes_cm':[lo*L/(nx-1),hi*L/(nx-1)],'marlpde_nodes_cm':[bounds[0]+(ml+.5)*dx,bounds[0]+(ml+1.5)*dx],'validation':checks,'rhythmite':diag(t,arrays['c_co'],arrays['u']),'rhythmite_bottom_vtu':diag(t,arrays['c_co_bottom'],arrays['u_bottom']),'rhythmite_bottom_ascii':diag(t,a[:,4],a[:,6]),'marlpde':diag(mt,mc)}
 rows.append(result)
 print(json.dumps(result),flush=True)
args.output_stem.with_suffix(".json").write_text(json.dumps({
    "measurement_depth_cm": 500,
    "method": "Linear interpolation; Marlpde L=500 linearly extrapolated from last two cell centres. Rhythmite uniform grid checked against rounded VTU coordinates; time from paired ASCII, not VTU plotting coordinates.",
    "rhythmite_tstar_years": args.rhythmite_tstar_years,
    "rhythmite_xstar_cm": args.rhythmite_xstar_cm,
    "rows": rows,
}, indent=2) + "\n")
keys = {"t_onset": "onset", "t_first_zero_U": "zero", "t_first_minimum_U": "minimum"}
times = {k: np.array([r["rhythmite"][v] if r["rhythmite"][v] is not None else np.nan for r in rows]) for k, v in keys.items()}
lengths = np.array([r["L"] for r in rows])
monsets = np.array([r["marlpde"]["onset"] if r["marlpde"]["onset"] is not None else np.nan for r in rows])
plot_times(lengths, times, lengths, monsets, output_stem=args.output_stem, show=False,
           measurement_label="Measurement depth: 500 cm")
