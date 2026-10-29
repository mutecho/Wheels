#include "Comparison.h"

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TPad.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>

namespace femtocompare {
namespace {

const std::array<const char*, 4> sampleNames{"A_SE", "A_ME", "B_SE", "B_ME"};
int blue, orange, sameColor, mixedColor;
unsigned long objectNumber = 0;

std::string Unique(const std::string& prefix) { return prefix + "_" + std::to_string(++objectNumber); }

void Text(double x, double y, const std::string& message, double size = 0.035, int color = 1) {
  TLatex t;
  t.SetNDC(); t.SetTextFont(42); t.SetTextSize(size); t.SetTextColor(color);
  t.DrawLatex(x, y, message.c_str());
}

std::unique_ptr<TH1D> RootHistogram(const Histogram& h, const std::string& name, bool normalized) {
  auto out = std::make_unique<TH1D>(name.c_str(), "", h.axis.bins(), h.axis.edges.data());
  out->SetDirectory(nullptr); out->Sumw2();
  out->GetXaxis()->SetTitle(h.axis.title.c_str());
  out->GetYaxis()->SetTitle(normalized ? "1/N dN/dx" : "Pairs / bin");
  for (int b = 0; b < h.axis.bins(); ++b) {
    const Point p = normalized ? NormalizedBin(h, b, true) : Point{h.count[b], h.variance[b], true};
    out->SetBinContent(b + 1, p.value);
    out->SetBinError(b + 1, std::sqrt(p.variance));
  }
  return out;
}

/// File ratios and the double ratio are dimensionless ratios of normalized
/// bin probabilities. Cross-file / SE-ME covariance is unavailable, so these
/// diagonal error bars are descriptive count-based approximations only.
Point RatioBin(const FourSpectra& f, int bin, int which) {
  const auto se = Divide(NormalizedBin(f.histogram[2], bin), NormalizedBin(f.histogram[0], bin));
  const auto me = Divide(NormalizedBin(f.histogram[3], bin), NormalizedBin(f.histogram[1], bin));
  if (which == 0) return se;
  if (which == 1) return me;
  return Divide(me, se);
}

TGraphErrors RatioGraph(const FourSpectra& f, int which, const std::string& name) {
  TGraphErrors result;
  result.SetName(name.c_str());
  for (int b = 0; b < f.histogram[0].axis.bins(); ++b) {
    const auto p = RatioBin(f, b, which);
    if (!p.valid) continue;
    const auto& edges = f.histogram[0].axis.edges;
    const int n = result.GetN();
    result.SetPoint(n, (edges[b] + edges[b + 1]) / 2, p.value);
    result.SetPointError(n, 0, std::sqrt(p.variance));
  }
  return result;
}

void AxisStyle(TH1& h, bool compact = false) {
  h.SetStats(false);
  h.GetXaxis()->SetTitleSize(compact ? 0.115 : 0.048);
  h.GetXaxis()->SetLabelSize(compact ? 0.095 : 0.042);
  h.GetXaxis()->SetTitleOffset(compact ? 1.02 : 1.12);
  h.GetYaxis()->SetTitleSize(compact ? 0.105 : 0.045);
  h.GetYaxis()->SetLabelSize(compact ? 0.085 : 0.04);
  h.GetYaxis()->SetTitleOffset(compact ? 0.60 : 1.45);
  h.GetYaxis()->SetNdivisions(compact ? 505 : 510);
}

void PadStyle(bool compact = false) {
  gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.035);
  gPad->SetTopMargin(compact ? 0.08 : 0.12); gPad->SetBottomMargin(compact ? 0.31 : 0.15);
  gPad->SetTicks(1, 1);
}

void DensityPair(const FourSpectra& f, int event, const Options& o) {
  PadStyle();
  auto a = RootHistogram(f.histogram[event], Unique("a_density"), true);
  auto b = RootHistogram(f.histogram[event + 2], Unique("b_density"), true);
  AxisStyle(*a); AxisStyle(*b);
  a->SetLineColor(blue); a->SetLineWidth(3);
  b->SetLineColor(orange); b->SetMarkerColor(orange); b->SetMarkerStyle(20); b->SetMarkerSize(0.7);
  a->SetMinimum(0); a->SetMaximum(1.34 * std::max({a->GetMaximum(), b->GetMaximum(), 1.e-10}));
  auto* drawnA = a->DrawCopy("HIST");
  auto* drawnB = b->DrawCopy("E1 SAME");
  TLegend legend(0.53, 0.68, 0.95, 0.86);
  legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42); legend.SetTextSize(0.036);
  legend.AddEntry(drawnA, o.labelA.c_str(), "l"); legend.AddEntry(drawnB, o.labelB.c_str(), "lep");
  legend.DrawClone();
  Text(0.16, 0.92, event == 0 ? "Same event" : "Mixed event", 0.047);
  std::ostringstream totals;
  totals << std::scientific << std::setprecision(3) << "N_{A}=" << f.histogram[event].total()
         << ", N_{B}=" << f.histogram[event + 2].total();
  Text(0.17, 0.80, totals.str(), 0.033);
}

void Ratios(const FourSpectra& f, bool doubleRatio, bool compact = false) {
  PadStyle(compact);
  auto frame = RootHistogram(f.histogram[0], Unique("ratio_frame"), false);
  frame->Reset(); AxisStyle(*frame, compact);
  frame->GetYaxis()->SetTitle(doubleRatio ? "(B/A)_{ME} / (B/A)_{SE}" : "B/A (unit area)");
  auto se = RatioGraph(f, doubleRatio ? 2 : 0, Unique("se_ratio"));
  auto me = RatioGraph(f, 1, Unique("me_ratio"));
  double low = 1, high = 1;
  for (auto* g : (doubleRatio ? std::vector<TGraphErrors*>{&se} : std::vector<TGraphErrors*>{&se, &me})) {
    for (int i = 0; i < g->GetN(); ++i) {
      low = std::min(low, g->GetPointY(i) - g->GetErrorY(i));
      high = std::max(high, g->GetPointY(i) + g->GetErrorY(i));
    }
  }
  const double span = std::max(0.08, high - low);
  frame->SetMinimum(std::max(0., low - 0.17 * span));
  frame->SetMaximum(high + 0.40 * span);
  frame->DrawCopy("AXIS");
  TLine unity(frame->GetXaxis()->GetXmin(), 1, frame->GetXaxis()->GetXmax(), 1);
  unity.SetLineStyle(2); unity.SetLineColor(kGray + 2); unity.DrawClone();
  se.SetMarkerStyle(20); se.SetMarkerSize(compact ? 0.55 : 0.75);
  se.SetMarkerColor(doubleRatio ? kBlack : sameColor); se.SetLineColor(doubleRatio ? kBlack : sameColor);
  auto* seDrawn = se.DrawClone("P SAME");
  if (!doubleRatio) {
    me.SetMarkerStyle(25); me.SetMarkerSize(compact ? 0.65 : 0.85);
    me.SetMarkerColor(mixedColor); me.SetLineColor(mixedColor);
    auto* meDrawn = me.DrawClone("P SAME");
    TLegend legend(0.18, compact ? 0.78 : 0.77, 0.72, compact ? 0.96 : 0.87);
    legend.SetNColumns(2); legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42);
    legend.SetTextSize(compact ? 0.09 : 0.035);
    legend.AddEntry(seDrawn, "Same B/A", "p"); legend.AddEntry(meDrawn, "Mixed B/A", "p"); legend.DrawClone();
  } else Text(0.17, 0.92, "Double ratio (unity = equal shape change)", 0.039);
  if (!compact) {
    const auto bins = f.histogram[0].axis.bins();
    const int missing = bins - se.GetN();
    if (missing) Text(0.17, 0.81, "Undefined bins omitted: " + std::to_string(missing), 0.031);
  }
  gPad->RedrawAxis();
}

void OverviewPanel(const FourSpectra& f) {
  auto* parent = gPad;
  auto* top = new TPad(Unique("top").c_str(), "", 0, 0.32, 1, 1);
  top->Draw(); top->cd();
  PadStyle(); gPad->SetBottomMargin(0.015); gPad->SetTopMargin(0.12);
  double maximum = 0;
  std::array<std::unique_ptr<TH1D>, 4> hist;
  for (int i = 0; i < 4; ++i) {
    hist[i] = RootHistogram(f.histogram[i], Unique("shape"), true);
    maximum = std::max(maximum, hist[i]->GetMaximum());
  }
  TLegend legend(0.20, 0.64, 0.93, 0.88);
  legend.SetNColumns(2); legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42); legend.SetTextSize(0.045);
  for (int i = 0; i < 4; ++i) {
    auto& h = hist[i]; AxisStyle(*h);
    h->GetXaxis()->SetLabelSize(0); h->GetXaxis()->SetTitleSize(0);
    h->GetYaxis()->SetTitleOffset(1.1); h->GetYaxis()->SetTitleSize(0.057); h->GetYaxis()->SetLabelSize(0.05);
    h->SetLineColor(i < 2 ? blue : orange); h->SetLineStyle(i % 2 == 0 ? 2 : 1); h->SetLineWidth(2);
    h->SetMinimum(0); h->SetMaximum((maximum > 0 ? maximum : 1.) * 1.47);
    auto* copy = h->DrawCopy(i == 0 ? "HIST" : "HIST SAME");
    const std::string label = std::string(i < 2 ? "A " : "B ") + (i % 2 == 0 ? "same" : "mixed");
    legend.AddEntry(copy, label.c_str(), "l");
  }
  legend.DrawClone(); Text(0.16, 0.92, f.histogram[0].axis.title, 0.056);
  if (maximum == 0) Text(0.21, 0.45, "No pairs in this selection", 0.060);
  parent->cd();
  auto* bottom = new TPad(Unique("bottom").c_str(), "", 0, 0, 1, 0.33);
  bottom->Draw(); bottom->cd(); Ratios(f, false, true);
  parent->cd();
}

FourSpectra GetSpectra(const Plan& p, const std::array<Sample, 4>& samples, int selection, int dimension) {
  FourSpectra f;
  f.name = p.selections[selection].name + "__" + p.axes[dimension].name;
  f.title = p.selections[selection].title + " | " + p.axes[dimension].name;
  for (int i = 0; i < 4; ++i) f.histogram[i] = samples[i].spectra[selection][dimension];
  return f;
}

TPad* SubPad(TCanvas& c, double x1, double y1, double x2, double y2) {
  c.cd();
  auto* pad = new TPad(Unique("panel").c_str(), "", x1, y1, x2, y2);
  pad->Draw(); pad->cd();
  return pad;
}

void Footer(TCanvas& canvas, const Options& o, const std::string& title) {
  canvas.cd();
  Text(0.025, 0.965, title, 0.028);
  Text(0.025, 0.026, o.labelA + "; " + o.labelB + ". Normal bins only. Errors: count approximation; pair/cross-file covariance unavailable.", 0.017);
}

/// Each plotted ratio also gets a complete-bin validity mask in the ROOT file;
/// undefined denominator bins are never encoded as valid 0 +/- 0 points.
void SaveSpectra(const FourSpectra& f, TDirectory& directory, std::ofstream& metrics, std::ofstream& bins) {
  directory.cd();
  for (int i = 0; i < 4; ++i) {
    RootHistogram(f.histogram[i], sampleNames[i], false)->Write();
    RootHistogram(f.histogram[i], std::string(sampleNames[i]) + "_density", true)->Write();
  }
  for (int r = 0; r < 3; ++r) {
    const std::string name = r == 0 ? "SE_B_over_A" : (r == 1 ? "ME_B_over_A" : "double_ME_over_SE");
    auto graph = RatioGraph(f, r, name); graph.Write();
    const auto& a = f.histogram[0].axis;
    TH1I mask((name + "_valid").c_str(), "1=valid, 0=undefined denominator", a.bins(), a.edges.data());
    mask.SetDirectory(nullptr);
    for (int b = 0; b < a.bins(); ++b) mask.SetBinContent(b + 1, RatioBin(f, b, r).valid);
    mask.Write();
  }
  metrics << f.name;
  for (const auto& h : f.histogram) metrics << ',' << h.total();
  metrics << ',' << TotalVariation(f.histogram[0], f.histogram[2]) << ',' << TotalVariation(f.histogram[1], f.histogram[3]) << '\n';
  for (int b = 0; b < f.histogram[0].axis.bins(); ++b) {
    bins << f.name << ',' << b + 1 << ',' << f.histogram[0].axis.edges[b] << ',' << f.histogram[0].axis.edges[b + 1];
    for (const auto& h : f.histogram) bins << ',' << h.count[b];
    for (int r = 0; r < 3; ++r) {
      const auto point = RatioBin(f, b, r);
      bins << ',' << point.valid << ',';
      if (point.valid) bins << point.value;
      bins << ',';
      if (point.valid) bins << std::sqrt(point.variance);
    }
    bins << '\n';
  }
}

std::unique_ptr<TH2D> RootMap(const Map2D& m, const std::string& name, bool normalized) {
  auto h = std::make_unique<TH2D>(name.c_str(), "", m.x.bins(), m.x.edges.data(), m.y.bins(), m.y.edges.data());
  h->SetDirectory(nullptr); h->Sumw2();
  h->GetXaxis()->SetTitle(m.x.title.c_str()); h->GetYaxis()->SetTitle(m.y.title.c_str());
  h->GetZaxis()->SetTitle(normalized ? "1/N d^{2}N/dxdy" : "Pairs / bin");
  const double total = std::accumulate(m.count.begin(), m.count.end(), 0.);
  const double totalV = std::accumulate(m.variance.begin(), m.variance.end(), 0.);
  for (int y = 0; y < m.y.bins(); ++y) for (int x = 0; x < m.x.bins(); ++x) {
    const int i = y * m.x.bins() + x;
    double value = m.count[i], variance = m.variance[i];
    if (normalized) {
      const double area = (m.x.edges[x + 1] - m.x.edges[x]) * (m.y.edges[y + 1] - m.y.edges[y]);
      const double p = total > 0 ? value / total : 0;
      variance = total > 0 ? std::max(0., ((1 - 2 * p) * variance + p * p * totalV) / (total * total * area * area)) : 0;
      value = p / area;
    }
    h->SetBinContent(x + 1, y + 1, value);
    h->SetBinError(x + 1, y + 1, std::sqrt(variance));
  }
  return h;
}

void DrawMap(TH2& h, const std::string& title, bool logz) {
  gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.19);
  gPad->SetTopMargin(0.15); gPad->SetBottomMargin(0.15);
  gPad->SetLogz(logz);
  h.SetStats(false);
  h.GetXaxis()->SetTitleSize(0.045); h.GetYaxis()->SetTitleSize(0.045); h.GetZaxis()->SetTitleSize(0.04);
  h.GetYaxis()->SetTitleOffset(1.5); h.GetZaxis()->SetTitleOffset(1.6);
  h.GetXaxis()->SetLabelSize(0.037); h.GetYaxis()->SetLabelSize(0.037); h.GetZaxis()->SetLabelSize(0.033);
  h.DrawCopy("COLZ");
  Text(0.15, 0.93, title, 0.04);
}

void MapsPage(TCanvas& canvas, const FourMaps& f, const Options& o, TDirectory& dir) {
  for (int event = 0; event < 2; ++event) {
    auto a = RootMap(f.histogram[event], Unique("map_a"), true);
    auto b = RootMap(f.histogram[event + 2], Unique("map_b"), true);
    const double maximum = std::max({1.e-10, a->GetMaximum(), b->GetMaximum()});
    a->SetMaximum(maximum * 1.02); b->SetMaximum(maximum * 1.02);
    a->SetMinimum(0); b->SetMinimum(0);
    dir.cd();
    RootMap(f.histogram[event], sampleNames[event], false)->Write();
    RootMap(f.histogram[event + 2], sampleNames[event + 2], false)->Write();
    const std::string ratioName = event ? "ME_B_over_A" : "SE_B_over_A";
    std::unique_ptr<TH2D> ratio(static_cast<TH2D*>(a->Clone(ratioName.c_str())));
    ratio->Reset(); ratio->GetZaxis()->SetTitle("B/A (unit area)");
    ratio->SetMinimum(); ratio->SetMaximum();
    std::unique_ptr<TH2D> mask(static_cast<TH2D*>(a->Clone((ratioName + "_valid").c_str())));
    mask->Reset(); mask->SetMinimum(0); mask->SetMaximum(1); mask->GetZaxis()->SetTitle("1=valid, 0=undefined");
    for (int y = 1; y <= a->GetNbinsY(); ++y) for (int x = 1; x <= a->GetNbinsX(); ++x) {
      const auto point = Divide({b->GetBinContent(x, y), std::pow(b->GetBinError(x, y), 2), true},
                                {a->GetBinContent(x, y), std::pow(a->GetBinError(x, y), 2), true});
      mask->SetBinContent(x, y, point.valid);
      if (point.valid) { ratio->SetBinContent(x, y, point.value); ratio->SetBinError(x, y, std::sqrt(point.variance)); }
    }
    ratio->Write(); mask->Write();
    const double y1 = event ? 0.08 : 0.51, y2 = event ? 0.50 : 0.93;
    const auto eventName = event ? "ME" : "SE";
    SubPad(canvas, 0.005, y1, 0.335, y2); DrawMap(*a, o.labelA + " " + eventName, false);
    SubPad(canvas, 0.335, y1, 0.665, y2); DrawMap(*b, o.labelB + " " + eventName, false);
    SubPad(canvas, 0.665, y1, 0.995, y2);
    // Fix the display range around unity so isolated low-count cells cannot
    // flatten the useful contrast. The complete, unclipped map was saved above.
    int clipped = 0;
    for (int y = 1; y <= mask->GetNbinsY(); ++y) for (int x = 1; x <= mask->GetNbinsX(); ++x) {
      if (mask->GetBinContent(x, y) &&
          (ratio->GetBinContent(x, y) < 0.8 || ratio->GetBinContent(x, y) > 1.2)) ++clipped;
    }
    ratio->SetMinimum(0.8); ratio->SetMaximum(1.2);
    DrawMap(*ratio, std::string(eventName) + " B/A; color 0.8-1.2 (" + std::to_string(clipped) + " clipped)", false);
    Text(0.16, 0.86, "Hatched: A=0; full values in comparison.root", 0.030);
    for (int y = 1; y <= mask->GetNbinsY(); ++y) for (int x = 1; x <= mask->GetNbinsX(); ++x) {
      if (mask->GetBinContent(x, y)) continue;
      TBox box(mask->GetXaxis()->GetBinLowEdge(x), mask->GetYaxis()->GetBinLowEdge(y),
               mask->GetXaxis()->GetBinUpEdge(x), mask->GetYaxis()->GetBinUpEdge(y));
      box.SetFillColor(kGray + 1); box.SetFillStyle(3004); box.SetLineColor(kGray + 1); box.DrawClone();
    }
  }
  Footer(canvas, o, f.title);
}

void DetailPage(TCanvas& canvas, const FourSpectra& f, const Options& o) {
  SubPad(canvas, 0.015, 0.51, 0.495, 0.94); DensityPair(f, 0, o);
  SubPad(canvas, 0.505, 0.51, 0.985, 0.94); DensityPair(f, 1, o);
  SubPad(canvas, 0.015, 0.075, 0.495, 0.50); Ratios(f, false);
  SubPad(canvas, 0.505, 0.075, 0.985, 0.50); Ratios(f, true);
  Footer(canvas, o, f.title);
}

void SummaryPage(TCanvas& canvas, const Options& o, const Plan& p, const std::array<Sample, 4>& samples, int selection) {
  const std::vector<int> dims = selection <= 1 ? std::vector<int>{0, 1, 2, 3, 4, 6} : std::vector<int>{0, 1, 2, 6};
  const int cols = dims.size() == 6 ? 3 : 2;
  for (std::size_t i = 0; i < dims.size(); ++i) {
    const int col = i % cols, row = i / cols;
    const double x1 = 0.01 + col * 0.98 / cols, x2 = 0.01 + (col + 1) * 0.98 / cols;
    const double y2 = 0.94 - row * 0.43, y1 = y2 - 0.425;
    SubPad(canvas, x1, y1, x2, y2);
    OverviewPanel(GetSpectra(p, samples, selection, dims[i]));
  }
  Footer(canvas, o, p.selections[selection].title + " | normalized spectra and B/A");
}

Point RawYieldBin(const FourSpectra& f, int bin, int which) {
  const auto raw = [&](int i) { return Point{f.histogram[i].count[bin], f.histogram[i].variance[bin], true}; };
  const auto a = Divide(raw(1), raw(0)), b = Divide(raw(3), raw(2));
  return which == 0 ? a : (which == 1 ? b : Divide(b, a));
}

void RawYieldPanel(const FourSpectra& f, bool relative, TDirectory& dir, std::ofstream& csv) {
  PadStyle();
  auto frame = RootHistogram(f.histogram[0], Unique("yield_frame"), false);
  frame->Reset(); AxisStyle(*frame);
  frame->GetYaxis()->SetTitle(relative ? "(ME/SE)_{B} / (ME/SE)_{A}" : "Raw pair yield ME / SE");
  std::array<TGraphErrors, 2> graph;
  const int curves = relative ? 1 : 2;
  double low = relative ? 1 : 1.e100, high = relative ? 1 : 0;
  for (int k = 0; k < curves; ++k) {
    const int which = relative ? 2 : k;
    const std::string name = relative ? "raw_yield_B_over_A" : (k == 0 ? "A_ME_over_SE" : "B_ME_over_SE");
    graph[k].SetName(name.c_str());
    const auto& axis = f.histogram[0].axis;
    TH1I mask((name + "_valid").c_str(), "", axis.bins(), axis.edges.data());
    mask.SetDirectory(nullptr);
    for (int b = 0; b < axis.bins(); ++b) {
      const auto point = RawYieldBin(f, b, which);
      mask.SetBinContent(b + 1, point.valid);
      csv << f.name << ',' << name << ',' << axis.edges[b] << ',' << axis.edges[b + 1] << ',' << point.valid << ',';
      if (point.valid) csv << point.value;
      csv << ',';
      if (point.valid) csv << std::sqrt(point.variance);
      csv << '\n';
      if (!point.valid) continue;
      int n = graph[k].GetN();
      graph[k].SetPoint(n, (axis.edges[b] + axis.edges[b + 1]) / 2, point.value);
      graph[k].SetPointError(n, 0, std::sqrt(point.variance));
      low = std::min(low, point.value - std::sqrt(point.variance));
      high = std::max(high, point.value + std::sqrt(point.variance));
    }
    dir.cd(); graph[k].Write(); mask.Write();
  }
  const double span = std::max(0.08, high - low);
  frame->SetMinimum(std::max(0., low - span * 0.15)); frame->SetMaximum(high + span * 0.35);
  frame->DrawCopy("AXIS");
  if (relative) {
    TLine line(frame->GetXaxis()->GetXmin(), 1, frame->GetXaxis()->GetXmax(), 1);
    line.SetLineStyle(2); line.SetLineColor(kGray + 2); line.DrawClone();
  }
  TLegend legend(0.18, 0.78, 0.63, 0.86);
  legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42); legend.SetTextSize(0.035); legend.SetNColumns(2);
  for (int k = 0; k < curves; ++k) {
    graph[k].SetLineColor(relative ? kBlack : (k == 0 ? blue : orange));
    graph[k].SetMarkerColor(graph[k].GetLineColor()); graph[k].SetMarkerStyle(k == 0 ? 20 : 25); graph[k].SetMarkerSize(0.8);
    auto* drawn = graph[k].DrawClone("P SAME");
    legend.AddEntry(drawn, relative ? "Relative yield B/A" : (k == 0 ? "A" : "B"), "p");
  }
  legend.DrawClone();
  Text(0.16, 0.92, relative ? "Relative mixing yield (no area normalization)" : "Mixing yield (no area normalization)", 0.038);
}

void YieldSummary(TCanvas& canvas, const Options& o, const Plan& p, const std::array<Sample, 4>& samples,
                  TDirectory& root, std::ofstream& csv) {
  const auto cent = GetSpectra(p, samples, 0, 4), mt = GetSpectra(p, samples, 0, 3);
  auto* centDir = root.mkdir("centrality"); auto* mtDir = root.mkdir("mt");
  SubPad(canvas, 0.015, 0.51, 0.495, 0.94); RawYieldPanel(cent, false, *centDir, csv);
  SubPad(canvas, 0.505, 0.51, 0.985, 0.94); RawYieldPanel(cent, true, *centDir, csv);
  SubPad(canvas, 0.015, 0.075, 0.495, 0.50); RawYieldPanel(mt, false, *mtDir, csv);
  SubPad(canvas, 0.505, 0.075, 0.985, 0.50);
  Text(0.08, 0.88, "Counts and integrated ME/SE", 0.050);
  std::ostringstream a, b, all, normal;
  a << std::scientific << std::setprecision(4) << "Normal bins A: SE=" << samples[0].normalWeight << ", ME=" << samples[1].normalWeight;
  b << std::scientific << std::setprecision(4) << "Normal bins B: SE=" << samples[2].normalWeight << ", ME=" << samples[3].normalWeight;
  const double ra = samples[1].normalWeight / samples[0].normalWeight, rb = samples[3].normalWeight / samples[2].normalWeight;
  normal << std::fixed << std::setprecision(5) << "Normal-bin ME/SE: A=" << ra << ", B=" << rb << ", B/A=" << rb / ra;
  const double aa = samples[1].allWeight / samples[0].allWeight, ab = samples[3].allWeight / samples[2].allWeight;
  all << std::fixed << std::setprecision(5) << "All-bin ME/SE: A=" << aa << ", B=" << ab << ", B/A=" << ab / aa;
  Text(0.08, 0.72, a.str(), 0.037); Text(0.08, 0.60, b.str(), 0.037);
  Text(0.08, 0.46, normal.str(), 0.037); Text(0.08, 0.34, all.str(), 0.037);
  Text(0.08, 0.18, "All-bin totals include flows for accounting only.", 0.036);
  Text(0.08, 0.08, "These yield ratios are not correlation functions.", 0.036);
  Footer(canvas, o, "Raw mixed / same pair yield | native bins, no area normalization");
}

}  // namespace

void WritePlots(const Options& o, const Plan& p, const std::array<Sample, 4>& samples,
                const std::vector<FourSpectra>& qa, const std::vector<FourMaps>& qaMaps) {
  blue = TColor::GetColor("#234E75"); orange = TColor::GetColor("#D07325");
  sameColor = TColor::GetColor("#237A78"); mixedColor = TColor::GetColor("#9D3E6D");
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetPalette(kViridis);
  gStyle->SetCanvasColor(0); gStyle->SetPadColor(0); gStyle->SetFrameLineColor(kGray + 2);
  TH1::AddDirectory(false);
  std::filesystem::create_directories(o.output / "png");
  TFile root((o.output / "comparison.root").c_str(), "RECREATE");
  if (root.IsZombie()) throw std::runtime_error("Cannot create output ROOT file");
  TNamed policy("comparison_policy", "All seven axes normal bins; no rebinning; density has unit area; R=shape_B/shape_A; D=R_ME/R_SE; errors approximate; validity masks accompany ratios.");
  policy.Write();
  std::ofstream metrics(o.output / "metrics.csv"), bins(o.output / "bin_values.csv"), index(o.output / "plot_index.csv");
  metrics << "distribution,N_A_SE,N_A_ME,N_B_SE,N_B_ME,TV_SE,TV_ME\n" << std::setprecision(16);
  bins << "distribution,bin,low,high,A_SE,A_ME,B_SE,B_ME,SE_valid,SE_BA,SE_error,ME_valid,ME_BA,ME_error,double_valid,double_ratio,double_error\n" << std::setprecision(16);
  index << "pdf_page,png,description\n";
  const auto pdf = (o.output / "comparison.pdf").string();
  TCanvas book("comparison_book", "", 1800, 1100);
  book.Print((pdf + "[").c_str());
  int pages = 0;
  auto save = [&](TCanvas& canvas, const std::string& name, const std::string& title) {
    const auto image = (o.output / "png" / (name + ".png")).string();
    canvas.SaveAs(image.c_str());
    canvas.Print(pdf.c_str(), ("Title:" + name).c_str());
    index << ++pages << ",png/" << name << ".png,\"" << title << "\"\n";
  };
  // Serialize every native-bin slice, including those not printed as a
  // redundant projection in the compact sliced overview pages.
  auto* spectraRoot = root.mkdir("spectra");
  for (std::size_t s = 0; s < p.selections.size(); ++s) {
    auto* selectionDir = spectraRoot->mkdir(p.selections[s].name.c_str());
    for (int d = 0; d < dimensions; ++d) {
      const auto f = GetSpectra(p, samples, s, d);
      SaveSpectra(f, *selectionDir->mkdir(p.axes[d].name.c_str()), metrics, bins);
    }
  }
  {
    TCanvas canvas(Unique("overview").c_str(), "", 1800, 1100);
    SummaryPage(canvas, o, p, samples, 0); save(canvas, "00_overview", "All-normal-bin overview");
  }
  {
    std::ofstream yields(o.output / "pair_yields.csv");
    yields << "distribution,ratio,low,high,valid,value,error\n" << std::setprecision(16);
    TCanvas canvas(Unique("yield_summary").c_str(), "", 1600, 1050);
    YieldSummary(canvas, o, p, samples, *root.mkdir("relative_mixing_yield"), yields);
    save(canvas, "00_pair_yields", "Raw mixed/same yield, centrality and mT");
  }
  for (int d = 0; d < dimensions; ++d) {
    TCanvas canvas(Unique("detail").c_str(), "", 1600, 1050);
    const auto f = GetSpectra(p, samples, 0, d);
    DetailPage(canvas, f, o); save(canvas, "01_" + p.axes[d].name, f.title);
  }
  for (std::size_t s = 1; s < p.selections.size(); ++s) {
    TCanvas canvas(Unique("slice").c_str(), "", p.selections[s].lowMomentum ? 1800 : 1500, 1100);
    SummaryPage(canvas, o, p, samples, s); save(canvas, "10_" + p.selections[s].name, p.selections[s].title);
  }
  auto* mapRoot = root.mkdir("maps");
  for (std::size_t k = 0; k < p.mapAxes.size(); ++k) {
    const auto pair = p.mapAxes[k];
    FourMaps maps{p.axes[pair[0]].name + "_vs_" + p.axes[pair[1]].name, "", {}};
    maps.title = "All-normal-bin 2D comparison | " + maps.name;
    for (int i = 0; i < 4; ++i) maps.histogram[i] = samples[i].maps[k];
    TCanvas canvas(Unique("maps").c_str(), "", 1800, 1100);
    MapsPage(canvas, maps, o, *mapRoot->mkdir(maps.name.c_str()));
    save(canvas, "20_" + maps.name, maps.title);
  }
  auto* qaRoot = root.mkdir("close_pair_qa");
  for (const auto& f : qa) {
    SaveSpectra(f, *qaRoot->mkdir(f.name.c_str()), metrics, bins);
    TCanvas canvas(Unique("qa").c_str(), "", 1600, 1050);
    DetailPage(canvas, f, o); save(canvas, "30_" + f.name, f.title);
  }
  for (const auto& maps : qaMaps) {
    TCanvas canvas(Unique("qa_map").c_str(), "", 1800, 1100);
    MapsPage(canvas, maps, o, *qaRoot->mkdir((maps.name + "_2d").c_str()));
    save(canvas, "31_" + maps.name, maps.title);
  }
  book.Print((pdf + "]").c_str());
  root.Close();
  std::ofstream report(o.output / "summary.md");
  report << "# Same / mixed comparison\n\n"
         << "- A: " << o.fileA << "\n- B: " << o.fileB << "\n"
         << "- Task A: " << o.taskA << "\n- Task B: " << o.taskB << "\n"
         << "- " << pages << " pages in comparison.pdf; PNG names and page order in plot_index.csv.\n\n"
         << "All plotted sparse distributions require normal bins on all seven axes. "
         << "Inventory and per-axis flow weights are in inventory.csv and flows.csv. Flow categories overlap; do not add them. "
         << "Native bins are retained, and variable-width spectra are shown as densities with unit integral. "
         << "Each centrality / mT slice is normalized separately. R_SE and R_ME are B/A shape ratios. "
         << "D=R_ME/R_SE is a shape diagnostic, without high-q correlation-function normalization.\n\n"
         << "Raw integrated ME/SE in the all-normal-bin domain: A=" << samples[1].normalWeight / samples[0].normalWeight
         << ", B=" << samples[3].normalWeight / samples[2].normalWeight
         << ", relative B/A=" << (samples[3].normalWeight / samples[2].normalWeight) / (samples[1].normalWeight / samples[0].normalWeight)
         << ". See png/00_pair_yields.png and pair_yields.csv; these unnormalized yield ratios are not correlation functions.\n\n"
         << "Input Sumw2 is preserved. Without it, bin variances use the ROOT unweighted-count convention. "
         << "Normalization includes Cov(n_i,N)=Var(n_i). Pair reuse, possible overlap between files, and SE/ME covariance "
         << "cannot be recovered from these histograms; error bars are descriptive approximations, not significance tests. "
         << "Zero-denominator ratio bins are omitted and have explicit validity masks in comparison.root. "
         << "2D undefined bins are hatched. The displayed 2D ratio color range is 0.8-1.2, with the count of clipped bins printed; "
         << "all numerical ratios are saved without clipping. QA plots have their own native normal-bin domain and no 7D cuts.\n\n"
         << "The qn axis has " << p.axes[5].bins() << " bin(s). Additional tasks are listed in inputs.txt and are never pooled.\n\n"
         << "| Variable | TV same | TV mixed |\n|---|---:|---:|\n" << std::fixed << std::setprecision(5);
  for (int d : {0, 1, 2, 3, 4, 6}) {
    const auto f = GetSpectra(p, samples, 0, d);
    report << '|' << p.axes[d].name << '|' << TotalVariation(f.histogram[0], f.histogram[2]) << '|'
           << TotalVariation(f.histogram[1], f.histogram[3]) << "|\n";
  }
  report << "\nTV=0.5*sum(|p_A-p_B|), a descriptive shape distance in [0,1]. "
         << "It includes statistical fluctuations and depends on the preserved native binning.\n";
  std::cout << "Wrote " << pages << " pages, ROOT projections/ratio masks, and numerical CSVs.\n";
}

}  // namespace femtocompare
