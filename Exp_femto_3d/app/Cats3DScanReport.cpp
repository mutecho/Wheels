#include "Cats3DModel.h"
#include "exp_femto_3d/CfNormalization.h"
#include "exp_femto_3d/Workflow.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TLegend.h"
#include "TNamed.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

namespace {
namespace fs = std::filesystem;
using exp_femto_3d::Cats3DConfig;
using exp_femto_3d::Cats3DLambda;
using exp_femto_3d::Cats3DLambdaInput;
using exp_femto_3d::Cats3DModel;
constexpr double kWindow = 0.04;  // Same transverse-bin average as the fit workflow, in GeV/c.

std::string Number(double x, int precision = 5) {
  if (!std::isfinite(x)) return "unavailable";
  std::ostringstream s; s << std::setprecision(precision) << x; return s.str();
}
std::string Json(const std::string &s) {
  std::ostringstream o; o << '"';
  for (unsigned char c : s) {
    if (c == '"' || c == '\\') o << '\\' << c;
    else if (c == '\n') o << "\\n";
    else if (c == '\r') o << "\\r";
    else if (c == '\t') o << "\\t";
    else if (c < 32) o << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
    else o << c;
  }
  return o.str() + '"';
}
std::string Json(double x) { return std::isfinite(x) ? Number(x, 17) : "null"; }
std::string Csv(const std::string &s) {
  std::string out = "\""; for (char c : s) { if (c == '"') out += '"'; out += c; } return out + '"';
}
std::string Csv(double x) { return std::isfinite(x) ? Number(x, 17) : "NaN"; }

TTree &Tree(TDirectory &dir, const std::string &path) {
  auto *t = dynamic_cast<TTree *>(dir.Get(path.c_str()));
  if (!t) throw std::runtime_error("Missing tree: " + path);
  return *t;
}
std::string Named(TDirectory &dir, const std::string &path) {
  auto *n = dynamic_cast<TNamed *>(dir.Get(path.c_str()));
  if (!n) throw std::runtime_error("Missing metadata: " + path);
  return n->GetTitle();
}
double Leaf(TTree &t, const char *name) {
  auto *l = t.GetLeaf(name); if (!l) throw std::runtime_error(std::string("Missing branch: ") + name);
  return l->GetValue();
}
std::vector<std::string> Split(const std::string &s, char sep) {
  std::vector<std::string> parts; std::istringstream in(s); std::string p;
  while (std::getline(in, p, sep)) parts.push_back(p); return parts;
}
TDirectory &Directory(TDirectory &parent, const std::string &name) {
  auto *d = parent.GetDirectory(name.c_str()); if (!d) d = parent.mkdir(name.c_str());
  if (!d) throw std::runtime_error("Cannot create ROOT directory: " + name); return *d;
}

struct Slice {
  std::string id, group, cf_path;
  double cent_low, cent_high, mt_low, mt_high, phi;
  bool integrated;
  exp_femto_3d::LevyFitResult metadata;
};
struct Context {
  std::string contract, kernel, geometry, statistic, baseline;
  std::vector<std::string> names;
  std::vector<bool> fixed;
  Cats3DConfig cfg;
  Cats3DLambdaInput lambda;
  std::array<std::array<double, 3>, 3> axes{};  // nbins, lower edge, upper edge in GeV/c.
  double qmax = 0;
  bool full = false, use_baseline = false;
};
struct Point {
  Slice slice;
  std::string kind, status, reason, image, reconstruction_error;
  double coordinate = 0, objective = 0, delta = 0, edm = 0, wall_ms = 0;
  int index = 0, migrad = 0, hesse = 0, covariance = 0, ndf = 0, attempts = 0;
  bool valid = false, errors_valid = false, domain_valid = false, objective_valid = false;
  bool model_drawn = false, boundary = false;
  std::vector<double> values, errors;
  std::vector<int> lower, upper;
};
struct Report {
  std::string kind, stem;
  double coordinate;
  std::vector<std::size_t> points;
  std::vector<std::string> group_images;
  int valid = 0, errors = 0, pages = 0;
};

std::vector<Slice> ReadSlices(TFile &scan) {
  auto &t = Tree(scan, "meta/CATS3DSelectedSlices");
  TTreeReader reader(&t);
  TTreeReaderValue<std::string> id(reader,"slice_id"), group(reader,"group_id"), path(reader,"cf_object_path"),
      qn_label(reader,"qn_label"), mt_mode(reader,"mt_rebin_mode"), phi_mode(reader,"phi_rebin_mode");
  TTreeReaderValue<int> cent(reader,"centrality_index"), mt(reader,"mt_index"), qn(reader,"qn_index"),
      phi(reader,"phi_index"), mapped(reader,"build_uses_symmetric_phi_range"),
      mt_rebin(reader,"mt_rebin_enabled"), phi_rebin(reader,"phi_rebin_enabled"),
      qn_all(reader,"is_qn_integrated"), phi_all(reader,"is_phi_integrated");
  TTreeReaderValue<double> cl(reader,"cent_low"), ch(reader,"cent_high"), ml(reader,"mt_low"), mh(reader,"mt_high"),
      ql(reader,"qn_low"), qh(reader,"qn_high"), pl(reader,"raw_phi_low"), ph(reader,"raw_phi_high"),
      dl(reader,"display_phi_low"), dh(reader,"display_phi_high"), pc(reader,"display_phi_center");
  std::unique_ptr<TTreeReaderValue<std::string>> period_mode;
  if (t.GetBranch("phi_mode")) period_mode=std::make_unique<TTreeReaderValue<std::string>>(reader,"phi_mode");
  exp_femto_3d::CfNormalizationTreeReader normalization(t,reader);
  std::vector<Slice> out; std::set<std::string> ids;
  while (reader.Next()) {
    if (!ids.insert(*id).second || id->find('/') != std::string::npos || id->find("..") != std::string::npos)
      throw std::runtime_error("Duplicate or unsafe slice id: " + *id);
    Slice s{*id,*group,*path,*cl,*ch,*ml,*mh,*pc,*phi_all!=0,{}};
    auto &r=s.metadata; r.normalization=normalization.Read();
    r.phi_mode=period_mode ? **period_mode : "separate";
    r.slice_id=*id; r.group_id=*group; r.slice_directory="fits/"+*id;
    r.centrality_index=*cent; r.mt_index=*mt; r.qn_index=*qn; r.phi_index=*phi;
    r.cent_low=*cl; r.cent_high=*ch; r.mt_low=*ml; r.mt_high=*mh;
    r.qn_low=*ql; r.qn_high=*qh; r.qn_label=*qn_label;
    r.raw_phi_low=*pl; r.raw_phi_high=*ph; r.display_phi_low=*dl; r.display_phi_high=*dh;
    r.phi=*pc; r.fit_uses_symmetric_phi_range=*mapped!=0;
    r.mt_rebin_enabled=*mt_rebin!=0; r.mt_rebin_mode=*mt_mode;
    r.phi_rebin_enabled=*phi_rebin!=0; r.phi_rebin_mode=*phi_mode;
    r.is_qn_integrated=*qn_all!=0; r.is_phi_integrated=*phi_all!=0;
    out.push_back(std::move(s));
  }
  if(out.size()!=static_cast<std::size_t>(t.GetEntries()))throw std::runtime_error("Cannot read selected-slice catalog");
  return out;
}

/** Read the persisted contract rather than today's mutable TOML. Explicit native
 * cache loading is mandatory: report generation must never generate a new kernel.
 */
Context ReadContext(TFile &scan, const Slice &slice, bool require_kernel = true) {
  auto *dir = scan.GetDirectory(("inputs/" + slice.id).c_str());
  if (!dir) throw std::runtime_error("Missing scan input metadata: " + slice.id);
  Context c; c.contract = Named(*dir,"CATS3DContract"); c.kernel = Named(*dir,"CATS3DKernelFile");
  c.geometry = Named(*dir,"FitGeometry"); c.statistic = Named(*dir,"FitStatistic");
  c.baseline = Named(*dir,"Baseline"); c.full = c.geometry == "full";
  if (!c.full && c.geometry != "diag") throw std::runtime_error("Unsupported geometry: " + c.geometry);
  c.use_baseline = c.baseline == "1+b2*qnorm_GeV_squared";
  if (!c.use_baseline && c.baseline != "unity") throw std::runtime_error("Unsupported baseline: " + c.baseline);
  auto fields = Split(c.contract, '|');
  if (fields.size() != 15 || fields[0] != "native_cats3d_byq_v1" ||
      fields[1] != "q_GeV_to_MeV_no_half" || fields[2] != "k_table_qnorm_half")
    throw std::runtime_error("Unsupported persisted CATS3D contract: " + c.contract);
  for (int k = 0; k < 3; ++k) {
    auto a = Split(fields[3+k], ':'); if (a.size() != 3) throw std::runtime_error("Invalid persisted axis");
    for (int j=0;j<3;++j) c.axes[k][j] = std::stod(a[j]);
    double half_width = (c.axes[k][2]-c.axes[k][1])/(2*c.axes[k][0]);
    c.qmax = std::max({c.qmax, std::abs(c.axes[k][1]+half_width), std::abs(c.axes[k][2]-half_width)});
  }
  c.qmax += 1e-12;
  c.cfg.grid_bins = std::stoi(fields[6]); c.cfg.grid_min_fm = std::stod(fields[7]);
  c.cfg.grid_max_fm = std::stod(fields[8]); c.cfg.subdivisions = std::stoi(fields[9]);
  c.cfg.max_pw = std::stoi(fields[10]); c.cfg.with_coulomb = std::stoi(fields[11]) != 0;
  c.cfg.use_disk_cache = true; c.cfg.kernel_cache_file = c.kernel;
  if (require_kernel && !fs::is_regular_file(c.kernel)) throw std::runtime_error("Recorded kernel is missing: " + c.kernel);
  auto &parameters = Tree(*dir,"ParameterContract");
  std::string name; auto *pname = &name; parameters.SetBranchAddress("name", &pname);
  for (Long64_t i=0;i<parameters.GetEntries();++i) {
    parameters.GetEntry(i); c.names.push_back(name); c.fixed.push_back(Leaf(parameters,"fixed") != 0);
  }
  parameters.ResetBranchAddresses();
  auto &t = Tree(*dir,"LambdaInputMetadata");
  auto *pk=&c.lambda.kind, *pm=&c.lambda.meaning, *pf=&c.lambda.file;
  auto *ph=&c.lambda.histogram, *pu=&c.lambda.axis_unit;
  t.SetBranchAddress("kind",&pk); t.SetBranchAddress("meaning",&pm); t.SetBranchAddress("file",&pf);
  t.SetBranchAddress("histogram",&ph); t.SetBranchAddress("axis_unit",&pu); t.GetEntry(0);
  double input = Leaf(t,"input_value"), weak = Leaf(t,"weak_decay_fraction");
  if (std::isfinite(input)) c.lambda.value=input;
  if (std::isfinite(weak)) c.lambda.weak_decay_fraction=weak;
  t.ResetBranchAddresses(); return c;
}

std::string PointPath(const Point &p) { return "lambda_scan/"+p.kind+"/"+p.slice.group+"/"+p.slice.id; }

void ReadPoints(TFile &scan, const Slice &slice, const Context &context,
                const std::string &kind, std::vector<Point> &out) {
  std::string path="lambda_scan/"+kind+"/"+slice.group+"/"+slice.id;
  auto *t=dynamic_cast<TTree*>(scan.Get((path+"/ScanPoints").c_str())); if(!t) return;
  std::vector<double> values, errors; std::vector<int> lo,hi;
  auto *pv=&values,*pe=&errors; auto *pl=&lo,*ph=&hi;
  std::string status,reason; auto *ps=&status,*pr=&reason;
  t->SetBranchAddress("parameters",&pv); t->SetBranchAddress("errors",&pe);
  t->SetBranchAddress("at_lower_bound",&pl); t->SetBranchAddress("at_upper_bound",&ph);
  t->SetBranchAddress("status",&ps); t->SetBranchAddress("failure_reason",&pr);
  std::map<int,int> counts; auto &attempts=Tree(scan,path+"/AttemptPoints");
  for(Long64_t i=0;i<attempts.GetEntries();++i) {attempts.GetEntry(i);++counts[int(Leaf(attempts,"point_index"))];}
  for(Long64_t i=0;i<t->GetEntries();++i) {
    t->GetEntry(i); Point p; p.slice=slice; p.kind=kind;
    p.index=Leaf(*t,"point_index"); p.coordinate=Leaf(*t,"coordinate");
    p.objective=Leaf(*t,"objective"); p.delta=Leaf(*t,"delta_objective"); p.edm=Leaf(*t,"edm");
    p.ndf=Leaf(*t,"ndf"); p.migrad=Leaf(*t,"migrad_status"); p.hesse=Leaf(*t,"hesse_status");
    p.covariance=Leaf(*t,"covariance_quality"); p.wall_ms=Leaf(*t,"wall_ms");
    p.valid=Leaf(*t,"valid")!=0; p.errors_valid=Leaf(*t,"parameter_errors_valid")!=0;
    p.domain_valid=Leaf(*t,"model_domain_valid")!=0; p.objective_valid=Leaf(*t,"objective_valid")!=0;
    p.boundary=Leaf(*t,"at_scan_boundary")!=0; p.status=status; p.reason=reason;
    p.values=values; p.errors=errors; p.lower=lo; p.upper=hi; p.attempts=counts[p.index];
    if(p.valid && (p.values.size()!=context.names.size() || p.errors.size()!=context.names.size()))
      throw std::runtime_error("Persisted parameter layout mismatch: "+slice.id);
    // Setup failures may have no parameter vector. Preserve their page and use
    // unavailable values rather than dropping the point or inventing zeros.
    p.values.resize(context.names.size(),std::numeric_limits<double>::quiet_NaN());
    p.errors.resize(context.names.size(),std::numeric_limits<double>::quiet_NaN());
    out.push_back(std::move(p));
  }
  t->ResetBranchAddresses();
}

std::unique_ptr<TH3D> Grid(const Context &c) {
  auto h=std::make_unique<TH3D>("ComputedGrid","",int(c.axes[0][0]),c.axes[0][1],c.axes[0][2],
    int(c.axes[1][0]),c.axes[1][1],c.axes[1][2],int(c.axes[2][0]),c.axes[2][1],c.axes[2][2]);
  h->SetDirectory(nullptr); return h;
}

/** Average the same normal transverse bins for data and model. The visible
 * longitudinal axis may be wider for the data; model curves stop at computed bins.
 * ROOT's stored data errors are propagated by Projection*("e") then averaged.
 */
std::unique_ptr<TH1D> Projection(TH3D &data,const TH3D &grid,int axis,const std::string &name) {
  const TAxis *da[]={data.GetXaxis(),data.GetYaxis(),data.GetZaxis()};
  const TAxis *ga[]={grid.GetXaxis(),grid.GetYaxis(),grid.GetZaxis()};
  int first[]={1,1,1},last[]={1,1,1},count=1;
  for(int k=0;k<3;++k) if(k!=axis) {
    double low=std::max(-kWindow,ga[k]->GetXmin()),high=std::min(kWindow,ga[k]->GetXmax());
    first[k]=std::max(1,da[k]->FindFixBin(low+1e-9));
    last[k]=std::min(da[k]->GetNbins(),da[k]->FindFixBin(high-1e-9));
    if(last[k]<first[k]) throw std::runtime_error("No common projection bins");
    count*=last[k]-first[k]+1;
  }
  TH1D *h=axis==0 ? data.ProjectionX(name.c_str(),first[1],last[1],first[2],last[2],"e") :
          axis==1 ? data.ProjectionY(name.c_str(),first[0],last[0],first[2],last[2],"e") :
                    data.ProjectionZ(name.c_str(),first[0],last[0],first[1],last[1],"e");
  h->SetDirectory(nullptr); h->Scale(1.0/count); return std::unique_ptr<TH1D>(h);
}

std::string Label(const Slice &s) {
  return "cent "+Number(s.cent_low)+"-"+Number(s.cent_high)+"%, m_{T} "+Number(s.mt_low)+"-"+
    Number(s.mt_high)+" GeV/c^{2}, "+(s.integrated ? "#varphi integrated" : "#varphi="+Number(s.phi,3)+" rad");
}
std::string CoordinateLabel(const Report &r) {
  return r.kind=="constant" ? "fixed #lambda="+Number(r.coordinate) : "fixed input-table scale="+Number(r.coordinate);
}

void CopyRows(TFile &scan,const Point &p,TDirectory &dir) {
  TDirectory::TContext current(&dir);
  for(const char *name:{"ScanPoints","AttemptPoints"}) {
    auto &tree=Tree(scan,PointPath(p)+"/"+name);
    std::unique_ptr<TTree> copy(tree.CopyTree(("point_index=="+std::to_string(p.index)).c_str()));
    if(!copy) throw std::runtime_error("Cannot copy stored scan rows");
    copy->SetDirectory(nullptr); copy->Write(name,TObject::kOverwrite);
  }
}

void DrawPoint(Point &p,const Context &c,TH3D &data,TH3D &grid,TH3D *model,
               TDirectory &dir,const fs::path &output,const Report &report,std::size_t number) {
  TDirectory::TContext current(&dir);
  const char *axes[]={"X","Y","Z"};
  const char *titles[]={"q_{out} (GeV/c)","q_{side} (GeV/c)","q_{long} (GeV/c)"};
  TCanvas canvas("Canvas_FitDiagnostics","CATS3D saved scan point",1400,1000); canvas.Divide(2,2,.008,.008);
  std::array<std::unique_ptr<TH1D>,3> observed;
  std::array<std::unique_ptr<TGraph>,3> curves;
  std::array<std::unique_ptr<TLegend>,3> legends;
  for(int axis=0;axis<3;++axis) {
    canvas.cd(axis+1); gPad->SetLeftMargin(.14);gPad->SetBottomMargin(.14);gPad->SetTopMargin(.10);
    auto &h=observed[axis]; h=Projection(data,grid,axis,std::string("Data_Proj")+axes[axis]);
    h->SetTitle(("CATS3D Gaussian;"+std::string(titles[axis])+";C(#bf{q})").c_str());
    h->SetStats(false);h->SetLineColor(kBlue+2);h->SetMarkerColor(kBlue+2);h->SetMarkerStyle(20);h->SetMarkerSize(.65);
    h->GetXaxis()->SetRangeUser(c.axes[axis][1],c.axes[axis][2]);
    h->GetXaxis()->SetTitleSize(.045);h->GetYaxis()->SetTitleSize(.045);h->GetYaxis()->SetTitleOffset(1.3);
    h->GetXaxis()->SetLabelSize(.04);h->GetYaxis()->SetLabelSize(.04);
    double ymin=1,ymax=1;
    for(int b=1;b<=h->GetNbinsX();++b) if(std::abs(h->GetBinCenter(b))<=c.qmax && h->GetBinError(b)>0) {
      ymin=std::min(ymin,h->GetBinContent(b)-h->GetBinError(b));
      ymax=std::max(ymax,h->GetBinContent(b)+h->GetBinError(b));
    }
    if(model) {
      auto prediction=Projection(*model,grid,axis,std::string("Model_Proj")+axes[axis]);
      prediction->Write("",TObject::kOverwrite);
      curves[axis]=std::make_unique<TGraph>(prediction->GetNbinsX());
      curves[axis]->SetName((std::string("Fit_Proj")+axes[axis]).c_str());
      for(int b=1;b<=prediction->GetNbinsX();++b) {
        double y=prediction->GetBinContent(b);
        curves[axis]->SetPoint(b-1,prediction->GetBinCenter(b),y);ymin=std::min(ymin,y);ymax=std::max(ymax,y);
      }
      curves[axis]->SetLineColor(p.valid ? kOrange+7 : kGray+2);
      curves[axis]->SetLineStyle(p.valid ? 1 : 2); curves[axis]->SetLineWidth(3);
    }
    double span=std::max(.1,ymax-ymin);h->SetMinimum(ymin-.12*span);h->SetMaximum(ymax+.2*span);
    h->Draw("E1");if(curves[axis])curves[axis]->Draw("L SAME");
    legends[axis]=std::make_unique<TLegend>(.20,.77,.95,.89);
    legends[axis]->SetBorderSize(0);legends[axis]->SetFillStyle(0);legends[axis]->SetTextSize(.032);
    legends[axis]->AddEntry(h.get(),"Stored CF; transverse |q| < 0.04","lep");
    if(curves[axis])legends[axis]->AddEntry(curves[axis].get(),p.valid ? "Saved fitted model" : "Unconverged trial (diagnostic)","l");
    legends[axis]->Draw();h->Write("",TObject::kOverwrite);
    if(curves[axis])curves[axis]->Write("",TObject::kOverwrite);
  }
  canvas.cd(4);
  TPaveText text(.025,.025,.98,.98,"NDC");text.SetFillColor(0);text.SetBorderSize(0);
  text.SetTextFont(42);text.SetTextAlign(12);text.SetTextSize(.031);
  text.AddText(("CATS3D | Gaussian #alpha=2 | "+CoordinateLabel(report)).c_str());
  text.AddText(Label(p.slice).c_str());
  text.AddText(("Minimum: "+p.status+"; conditional errors: "+(p.errors_valid?"valid":"unavailable")).c_str());
  text.AddText(("MIGRAD="+std::to_string(p.migrad)+"; HESSE="+std::to_string(p.hesse)+
    "; covariance="+std::to_string(p.covariance)+"; EDM="+Number(p.edm,3)).c_str());
  text.AddText(("Stored objective="+Number(p.objective,7)+"; ndf="+std::to_string(p.ndf)+
    "; #Deltaobjective="+Number(p.delta,4)).c_str());
  for(std::size_t j=0;j<c.names.size();++j) if(!c.fixed[j]) {
    std::string v=c.names[j]+" = "+Number(p.values[j],6);
    v+=p.errors_valid && std::isfinite(p.errors[j]) ? " #pm "+Number(p.errors[j],3) : " (error unavailable)";
    if((j<p.lower.size()&&p.lower[j])||(j<p.upper.size()&&p.upper[j])) v+=" [BOUND]";
    text.AddText(v.c_str());
  }
  text.AddText(("Seed attempts="+std::to_string(p.attempts)+"; stored selected-attempt time="+Number(p.wall_ms/1000,4)+" s").c_str());
  if(!p.reason.empty())text.AddText(p.reason.c_str());
  if(!p.reconstruction_error.empty())text.AddText(("No model curve: "+p.reconstruction_error).c_str());
  text.AddText("Errors condition on fixed external #lambda and numerical kernel.");
  text.Draw();canvas.Modified();canvas.Update();
  p.image="figures/"+report.stem+"/point_"+std::to_string(number)+".png";
  canvas.Print((output/p.image).c_str());
  canvas.Write("",TObject::kOverwrite);
  if(model)model->Write("ReconstructedModelCF3D",TObject::kOverwrite);
  TNamed provenance("Reconstruction","Saved parameters; same CATS3D adapter and native kernel; no minimization; failed trials remain invalid");
  provenance.Write();
  canvas.Clear();  // Destroy pads before their non-owned histogram/graph primitives.
}

bool DrawGroup(Report &r,const std::string &group,const std::vector<Point> &points,
               const std::map<std::string,Context> &contexts,TDirectory &root,const fs::path &out) {
  std::vector<const Point*> members;
  for(auto index:r.points)if(points[index].slice.group==group && !points[index].slice.integrated)members.push_back(&points[index]);
  if(members.empty())return false;
  auto &dir=Directory(Directory(root,"groups"),group);TDirectory::TContext current(&dir);
  auto &context=contexts.at(members.front()->slice.id);
  std::vector<int> indices;for(std::size_t j=0;j<context.fixed.size();++j)if(!context.fixed[j])indices.push_back(j);
  TCanvas canvas("ParameterVsPhi","Saved CATS3D scan parameters vs phi",1400,1000);canvas.Divide(3,3,.01,.01);
  std::vector<std::unique_ptr<TGraphErrors>> measured;
  std::vector<std::unique_ptr<TGraph>> noerror;
  for(std::size_t k=0;k<indices.size();++k) {
    int j=indices[k];canvas.cd(k+1);gPad->SetLeftMargin(.18);gPad->SetBottomMargin(.15);
    auto good=std::make_unique<TGraphErrors>();auto missing=std::make_unique<TGraph>();
    double ymin=std::numeric_limits<double>::infinity(),ymax=-ymin;
    for(auto *p:members)if(p->valid && std::isfinite(p->values[j])) {
      double y=p->values[j];ymin=std::min(ymin,y);ymax=std::max(ymax,y);
      if(p->errors_valid && std::isfinite(p->errors[j])) {
        int n=good->GetN();good->SetPoint(n,p->slice.phi,y);good->SetPointError(n,0,p->errors[j]);
        ymin=std::min(ymin,y-p->errors[j]);ymax=std::max(ymax,y+p->errors[j]);
      }else missing->SetPoint(missing->GetN(),p->slice.phi,y);
    }
    if(!std::isfinite(ymin)){ymin=0;ymax=1;}
    double span=std::max(.05,ymax-ymin);
    auto *frame=gPad->DrawFrame(0,ymin-.15*span,3.141592653589793,ymax+.2*span);
    frame->SetTitle((context.names[j]+";#varphi (rad);"+context.names[j]).c_str());
    frame->GetXaxis()->SetTitleSize(.05);frame->GetYaxis()->SetTitleSize(.05);frame->GetYaxis()->SetTitleOffset(1.7);
    good->SetName((context.names[j]+"_with_errors").c_str());good->SetMarkerStyle(20);good->SetMarkerColor(kBlue+2);good->SetLineColor(kBlue+2);
    missing->SetName((context.names[j]+"_errors_unavailable").c_str());missing->SetMarkerStyle(24);missing->SetMarkerColor(kOrange+7);
    if(good->GetN())good->Draw("P SAME");if(missing->GetN())missing->Draw("P SAME");
    good->Write();missing->Write();measured.push_back(std::move(good));noerror.push_back(std::move(missing));
  }
  canvas.cd(9);TPaveText note(.04,.06,.96,.94,"NDC");note.SetBorderSize(0);note.SetFillColor(0);note.SetTextFont(42);note.SetTextSize(.047);note.SetTextAlign(12);
  note.AddText(CoordinateLabel(r).c_str());note.AddText(Label(members.front()->slice).c_str());
  note.AddText("Filled: valid minimum + conditional HESSE error");note.AddText("Open: valid minimum, error unavailable");
  note.AddText("Invalid minima omitted; phi-integrated result separate");note.AddText("No harmonic fit / automatic lambda selection");note.Draw();
  canvas.cd();canvas.Modified();canvas.Update();
  std::string image="figures/"+r.stem+"/group_"+std::to_string(r.group_images.size())+".png";
  canvas.Print((out/image).c_str());r.group_images.push_back(image);
  canvas.Write();canvas.Clear();
  return true;
}

void WriteManifest(const fs::path &out,const std::string &scan_path,const std::string &cf_path,
                   const std::vector<Report> &reports,const std::vector<Point> &points,
                   const std::map<std::string,Context> &contexts,double seconds,std::size_t folds) {
  std::ofstream f(out/"manifest.json");
  f<<"{\"complete\":true,\"scan_root\":"<<Json(scan_path)<<",\"input_cf_root\":"<<Json(cf_path)
   <<",\"elapsed_seconds\":"<<Json(seconds)<<",\"source_folds\":"<<folds<<",\"reports\":[";
  for(std::size_t i=0;i<reports.size();++i){const auto &r=reports[i];if(i)f<<',';
    f<<"{\"kind\":"<<Json(r.kind)<<",\"coordinate\":"<<Json(r.coordinate)<<",\"stem\":"<<Json(r.stem)
     <<",\"fit_report_root\":"<<Json("fit_report_"+r.stem+".root")
     <<",\"points\":"<<r.points.size()<<",\"valid\":"<<r.valid<<",\"errors_valid\":"<<r.errors<<",\"pages\":"<<r.pages
     <<",\"group_images\":[";
    for(std::size_t j=0;j<r.group_images.size();++j){if(j)f<<',';f<<Json(r.group_images[j]);}f<<"]}";}
  f<<"],\"points\":[";
  for(std::size_t i=0;i<points.size();++i){const auto &p=points[i];const auto &c=contexts.at(p.slice.id);if(i)f<<',';
    f<<"{\"id\":"<<Json(p.slice.id)<<",\"group\":"<<Json(p.slice.group)<<",\"kind\":"<<Json(p.kind)
     <<",\"coordinate\":"<<Json(p.coordinate)<<",\"phi_integrated\":"<<(p.slice.integrated?"true":"false")
     <<",\"phi\":"<<Json(p.slice.phi)<<",\"valid\":"<<(p.valid?"true":"false")<<",\"errors_valid\":"<<(p.errors_valid?"true":"false")
     <<",\"status\":"<<Json(p.status)<<",\"reason\":"<<Json(p.reason)<<",\"objective\":"<<Json(p.objective)
     <<",\"delta_objective\":"<<Json(p.delta)<<",\"ndf\":"<<p.ndf<<",\"edm\":"<<Json(p.edm)
     <<",\"migrad\":"<<p.migrad<<",\"hesse\":"<<p.hesse<<",\"covariance\":"<<p.covariance
     <<",\"attempts\":"<<p.attempts<<",\"model_drawn\":"<<(p.model_drawn?"true":"false")
     <<",\"reconstruction_error\":"<<Json(p.reconstruction_error)<<",\"image\":"<<Json(p.image)<<",\"parameters\":[";
    for(std::size_t j=0;j<c.names.size();++j){if(j)f<<',';f<<"{\"name\":"<<Json(c.names[j])<<",\"value\":"<<Json(p.values[j])
     <<",\"error\":"<<Json(p.errors[j])<<",\"fixed\":"<<(c.fixed[j]?"true":"false")<<",\"at_bound\":"
     <<(((j<p.lower.size()&&p.lower[j])||(j<p.upper.size()&&p.upper[j]))?"true":"false")<<'}';}
    f<<"]}";
  }
  f<<"]}\n";
}

/** Translate the persisted parameter order; lambda and alpha remain fixed inputs.
 * Failed minima and missing HESSE errors retain their saved values/statuses in
 * FitCatalog. The standard writer gates harmonic fits on complete usable groups.
 */
exp_femto_3d::LevyFitResult SummaryResult(const Point &p,const Context &c,const std::string &scan_path) {
  auto r=p.slice.metadata;
  r.fit_model=c.geometry; r.has_off_diagonal=c.full; r.correlation_backend="cats3d";
  r.uses_coulomb=c.cfg.with_coulomb; r.coulomb_mode="cats3d";
  r.uses_core_halo_lambda=false; r.uses_q2_baseline=c.use_baseline;
  r.uses_pml=c.statistic=="existing_raw_count_PML";
  r.lambda_mode="external_fixed"; r.lambda_is_external_fixed=true;
  r.lambda_input_kind=p.kind=="constant" ? "constant" : c.lambda.kind;
  r.lambda_input_meaning=p.kind=="constant" ? "final_lambda" : c.lambda.meaning;
  r.lambda_input_object_path=scan_path+":inputs/"+p.slice.id+"/LambdaInputMetadata";
  r.norm=p.values.at(0); r.norm_err=p.errors.at(0);
  r.lambda=p.kind=="constant" ? p.coordinate : std::numeric_limits<double>::quiet_NaN();
  r.lambda_err=p.kind=="constant" ? 0.0 : std::numeric_limits<double>::quiet_NaN();
  r.rout2=p.values.at(2); r.rout2_err=p.errors.at(2);
  r.rside2=p.values.at(3); r.rside2_err=p.errors.at(3);
  r.rlong2=p.values.at(4); r.rlong2_err=p.errors.at(4);
  if(c.full) {
    r.routside2=p.values.at(5); r.routside2_err=p.errors.at(5);
    r.routlong2=p.values.at(6); r.routlong2_err=p.errors.at(6);
    r.rsidelong2=p.values.at(7); r.rsidelong2_err=p.errors.at(7);
  }
  r.alpha=p.values.at(c.full?8:5); r.alpha_err=p.errors.at(c.full?8:5);
  r.baseline_q2=p.values.at(c.full?9:6); r.baseline_q2_err=p.errors.at(c.full?9:6);
  r.fit_statistic=p.objective; r.ndf=p.ndf; r.status=p.migrad; r.edm=p.edm; r.minuit_istat=p.covariance;
  r.point_estimate_valid=p.valid; r.parameter_errors_valid=p.errors_valid;
  r.fit_failure_reason=p.reason;
  r.at_parameter_boundary=std::any_of(p.lower.begin(),p.lower.end(),[](int v){return v!=0;}) ||
                          std::any_of(p.upper.begin(),p.upper.end(),[](int v){return v!=0;});
  return r;
}

void WriteStandardReport(TFile &scan,const fs::path &input,const fs::path &output,const Report &r,
                         const std::vector<Point> &points,const std::map<std::string,Context> &contexts) {
  std::vector<exp_femto_3d::LevyFitResult> results;
  for(auto index:r.points)results.push_back(SummaryResult(points[index],contexts.at(points[index].slice.id),input.string()));
  const auto path=output/("fit_report_"+r.stem+".root");
  exp_femto_3d::WriteCats3DFitReportRoot(path.string(),results);
  TFile file(path.c_str(),"UPDATE");if(file.IsZombie())throw std::runtime_error("Cannot reopen standard report");
  TNamed("InputScanRoot",input.c_str()).Write();TNamed("InputCFRoot",Named(scan,"InputCFRoot").c_str()).Write();
  TNamed("ReportMode","saved_scan_parameters; standard_phi_harmonics; no_slice_minimization_or_kernel").Write();
  TNamed("ScanCoordinate",CoordinateLabel(r).c_str()).Write();
  for(auto index:r.points) {
    const auto &p=points[index];auto &dir=Directory(Directory(file,"slices"),p.slice.id);CopyRows(scan,p,dir);
    TDirectory::TContext current(&dir);
    TNamed("CATS3DContract",contexts.at(p.slice.id).contract.c_str()).Write();
    const auto source="inputs/"+p.slice.id+"/";
    auto &meta=Tree(scan,source+"LambdaInputMetadata");
    std::unique_ptr<TTree> copy(meta.CloneTree(-1,"fast"));copy->SetDirectory(nullptr);copy->Write();
    // This is the original table, before the recorded coordinate's scale.
    if(auto *table=dynamic_cast<TH1*>(scan.Get((source+"LambdaInput").c_str())))table->Write("LambdaInput");
  }
  file.cd();TNamed("ReportStatus","complete").Write();file.Close();
  std::cout<<"[report] standard fit report: "<<path<<std::endl;
}

void Run(const fs::path &input,const fs::path &output,bool summary_only) {
  const auto start=std::chrono::steady_clock::now();
  if(!summary_only && !Cats3DModel::Available())throw std::runtime_error("Report reconstruction requires a CATS3D-enabled build");
  if(fs::exists(output))throw std::runtime_error("Output directory already exists; choose a new directory: "+output.string());
  TFile scan(input.c_str(),"READ");if(scan.IsZombie())throw std::runtime_error("Cannot open scan ROOT");
  if(Named(scan,"CATS3DRunStatus")!="complete")throw std::runtime_error("Scan is not marked complete");
  std::string cf_path=Named(scan,"InputCFRoot");
  std::unique_ptr<TFile> cf;
  if(!summary_only){cf=std::make_unique<TFile>(cf_path.c_str(),"READ");
    if(cf->IsZombie())throw std::runtime_error("Cannot open recorded CF ROOT: "+cf_path);}
  auto slices=ReadSlices(scan);std::map<std::string,Context> contexts;std::vector<Point> points;
  for(const auto &s:slices){auto c=ReadContext(scan,s,!summary_only);
    for(const std::string kind:{"constant","scale"})ReadPoints(scan,s,c,kind,points);
    contexts.emplace(s.id,std::move(c));}
  if(points.empty())throw std::runtime_error("No stored scan points");
  std::map<std::pair<std::string,double>,Report> grouped;
  for(std::size_t i=0;i<points.size();++i){auto &r=grouped[{points[i].kind,points[i].coordinate}];
    r.kind=points[i].kind;r.coordinate=points[i].coordinate;r.points.push_back(i);r.valid+=points[i].valid;r.errors+=points[i].valid&&points[i].errors_valid;}
  fs::create_directories(output);std::ofstream(output/"report_status.txt")<<"incomplete\n";
  // Snapshot table inputs from the scan itself, so changed or moved original
  // calibration files cannot silently alter reconstructed scale-scan models.
  if(!summary_only && std::any_of(contexts.begin(),contexts.end(),[](const auto &x){return x.second.lambda.kind=="histogram";})) {
    TFile tables((output/"lambda_inputs.root").c_str(),"CREATE");
    for(auto &[id,c]:contexts)if(c.lambda.kind=="histogram") {
      auto *h=dynamic_cast<TH1*>(scan.Get(("inputs/"+id+"/LambdaInput").c_str()));
      if(!h)throw std::runtime_error("Persisted lambda input histogram missing: "+id);
      tables.cd();h->Write(id.c_str());c.lambda.file=(output/"lambda_inputs.root").string();c.lambda.histogram=id;
    }
    tables.Close();
  }
  std::vector<Report> reports;std::unique_ptr<Cats3DModel> model;std::string active_contract;
  std::size_t folds=0,processed=0;
  for(auto &[key,r]:grouped) {
    std::string token=Number(r.coordinate,9);for(char &ch:token)if(ch=='.')ch='p';else if(ch=='-')ch='m';else if(ch=='+')ch='_';
    r.stem="cats3d_"+r.kind+"_"+std::to_string(reports.size())+"_"+token;
    if(r.points.size()!=slices.size())throw std::runtime_error("Scan coordinate has missing slices; refusing a reduced group summary");
    WriteStandardReport(scan,input,output,r,points,contexts);
    if(summary_only){reports.push_back(r);continue;}
    fs::create_directories(output/"figures"/r.stem);
    // The Python presentation step assembles these exact lossless ROOT images
    // into consistently sized PDF pages, independently of ROOT's paper backend.
    r.pages=1;
    TFile root((output/(r.stem+".root")).c_str(),"CREATE");if(root.IsZombie())throw std::runtime_error("Cannot create report ROOT");
    TNamed("InputScanRoot",input.c_str()).Write();TNamed("InputCFRoot",cf_path.c_str()).Write();
    TNamed("ReportMode","offline_saved_scan_parameters_no_minimization").Write();
    TNamed("ReportStatus","incomplete").Write();
    std::set<std::string> groups;for(auto index:r.points)groups.insert(points[index].slice.group);
    for(const auto &group:groups)r.pages+=DrawGroup(r,group,points,contexts,root,output);
    for(auto index:r.points) {
      auto &p=points[index];auto &c=contexts.at(p.slice.id);
      std::unique_ptr<TH3D> data(dynamic_cast<TH3D*>(cf->Get(p.slice.cf_path.c_str())));
      if(!data)throw std::runtime_error("Missing stored CF: "+p.slice.cf_path);data->SetDirectory(nullptr);
      auto grid=Grid(c);std::unique_ptr<TH3D> prediction;
      if(!model || c.contract!=active_contract) {
        if(model)folds+=model->FoldCount();model.reset();
        auto weights=std::make_shared<Cats3DLambda>(c.lambda);
        Cats3DModel::ValidateInputGrid(*data,c.qmax,c.cfg,*weights);
        model=std::make_unique<Cats3DModel>(*data,c.qmax,c.cfg);
        if(model->Contract()!=c.contract)throw std::runtime_error("Current adapter does not reproduce the persisted kernel contract");
        active_contract=c.contract;
        std::cout<<"[report] kernel loaded; resident kernel GiB="<<model->KernelBytes()/double(1ULL<<30)
          <<"; estimated process peak GiB="<<model->EstimatedPeakGiB()<<std::endl;
      }
      model->SetLambda(std::make_shared<Cats3DLambda>(c.lambda));
      if(p.domain_valid && p.objective_valid && std::all_of(p.values.begin(),p.values.end(),[](double x){return std::isfinite(x);})) {
        try{
          model->SetScan(p.kind=="constant" ? std::optional<double>(p.coordinate) : std::nullopt,p.kind=="scale"?p.coordinate:1.0);
          prediction=model->FinalHistogram(p.values.data(),c.full,c.use_baseline);p.model_drawn=true;
        }
        catch(const exp_femto_3d::Cats3DDomainError &error){p.reconstruction_error=error.what();}
      } else p.reconstruction_error="Saved model/objective domain is invalid";
      if(p.valid && !prediction)throw std::runtime_error("Cannot reconstruct a valid saved minimum: "+p.slice.id+" "+p.reconstruction_error);
      auto &dir=Directory(Directory(root,"slices"),p.slice.id);CopyRows(scan,p,dir);
      {TDirectory::TContext current(&dir);TNamed("CATS3DContract",c.contract.c_str()).Write();TNamed("CATS3DKernelFile",c.kernel.c_str()).Write();}
      DrawPoint(p,c,*data,*grid,prediction.get(),dir,output,r,index);++r.pages;
      if(++processed%10==0 || processed==points.size())std::cout<<"[report] rendered "<<processed<<"/"<<points.size()<<std::endl;
    }
    root.cd();TNamed("ReportStatus","complete").Write("",TObject::kOverwrite);root.Close();
    reports.push_back(r);
  }
  if(model)folds+=model->FoldCount();model.reset();
  std::ofstream csv(output/"scan_points.csv");
  csv<<"kind,coordinate,group_id,slice_id,valid,parameter_errors_valid,status,objective,ndf,delta_objective,edm,migrad_status,hesse_status,covariance_quality,attempts,parameter,value,error,fixed,at_bound\n";
  for(const auto &p:points){const auto &c=contexts.at(p.slice.id);for(std::size_t j=0;j<c.names.size();++j)
    csv<<Csv(p.kind)<<','<<Csv(p.coordinate)<<','<<Csv(p.slice.group)<<','<<Csv(p.slice.id)<<','<<p.valid<<','<<p.errors_valid<<','<<Csv(p.status)<<','
    <<Csv(p.objective)<<','<<p.ndf<<','<<Csv(p.delta)<<','<<Csv(p.edm)<<','<<p.migrad<<','<<p.hesse<<','<<p.covariance<<','<<p.attempts<<','
    <<Csv(c.names[j])<<','<<Csv(p.values[j])<<','<<Csv(p.errors[j])<<','<<c.fixed[j]<<','
    <<((j<p.lower.size()&&p.lower[j])||(j<p.upper.size()&&p.upper[j]))<<'\n';}
  double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  WriteManifest(output,input.string(),cf_path,reports,points,contexts,seconds,folds);
  std::ofstream(output/"report_status.txt")<<"complete\n";
  std::cout<<"[report] complete; points="<<points.size()<<" folds="<<folds<<" seconds="<<seconds<<" output="<<output<<std::endl;
}
}  // namespace

int main(int argc,char **argv) {
  try {
    fs::path scan,output;bool summary_only=false;
    for(int i=1;i<argc;++i) {
      std::string arg=argv[i];
      if(arg=="--help"||arg=="-h") {
        std::cout<<"Usage: cats3d_scan_report --scan-root FILE --output-directory NEW_DIRECTORY [--summary-only]\n"
          <<"Export ROOT/PNG and a manifest using the recorded native kernel; no minimization.\n"
          <<"--summary-only: standard fit_report ROOTs only, without CF/kernel loading or source folding.\n"
          <<"Use scripts/run_cats3d_scan_report.sh for the PDF and HTML presentation step.\n";return 0;
      }
      if(arg=="--summary-only"){summary_only=true;continue;}
      if((arg=="--scan-root"||arg=="--output-directory")&&i+1<argc) {
        if(arg=="--scan-root")scan=fs::absolute(argv[++i]);else output=fs::absolute(argv[++i]);
      }else throw std::runtime_error("Unknown option or missing value: "+arg);
    }
    if(scan.empty()||output.empty())throw std::runtime_error("Both --scan-root and --output-directory are required");
    gROOT->SetBatch(kTRUE);TH1::AddDirectory(kFALSE);gStyle->SetOptStat(0);gStyle->SetOptTitle(1);gStyle->SetTextFont(42);
    Run(scan,output,summary_only);return 0;
  }catch(const std::exception &error){std::cerr<<"[error] "<<error.what()<<std::endl;return 1;}
}
