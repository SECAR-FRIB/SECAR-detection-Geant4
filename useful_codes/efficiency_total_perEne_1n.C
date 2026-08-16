#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <algorithm>

// pack (eventID, trackID) into one key
static long long Key32(int eventID, int trackID) {
  return ( (long long)eventID << 32 ) | (unsigned int)trackID;
}

// old-style threshold logic
static inline bool PassThreshold(double edep, int Zdep, int Adep, double thr, double heavy_thr)
{
  if (Zdep == 1 && Adep == 1) return (edep > thr);
  if (Zdep > 1 && Adep > 1)   return (edep > heavy_thr);
  return false;
}

struct GraphPair {
  TGraphErrors* pts = nullptr;   // points with stat errors
  TGraphErrors* band = nullptr;  // threshold band
};

static GraphPair MakePair() { return { new TGraphErrors(), new TGraphErrors() }; }

// Use a positive requested limit verbatim. Otherwise derive a padded limit
// from both the statistical-error points and threshold-uncertainty band.
static double ResolveYMax(const GraphPair& graph, double requested)
{
  if (requested > 0.0) return requested;

  double plottedMax = 0.0;
  for (const TGraphErrors* item : {graph.pts, graph.band}) {
    if (!item) continue;
    for (int i = 0; i < item->GetN(); ++i) {
      double x = 0.0, y = 0.0;
      item->GetPoint(i, x, y);
      plottedMax = std::max(plottedMax, y + item->GetErrorY(i));
    }
  }

  // Leave headroom above the highest point. Keep an informative range for
  // empty/all-zero graphs instead of producing a degenerate axis.
  return plottedMax > 0.0 ? 1.15 * plottedMax : 1.0;
}

static void FillEventGraph(
  const std::vector<double>& x,
  const std::vector<int>& denomEvents,
  const std::vector<int>& num_nom,
  const std::vector<int>& num_max,
  const std::vector<int>& num_min,
  const std::vector<double>& xerrVec,
  GraphPair& out
){
  for (size_t b=0; b<x.size(); b++) {
    const double denom = (double)denomEvents[b];
    if (denom <= 0) continue;

    const double num = (double)num_nom[b];
    const double eff = 100.0*(num/denom);

    const double eff_max = 100.0*((double)num_max[b]/denom);
    const double eff_min = 100.0*((double)num_min[b]/denom);
    const double thrErr  = std::fabs(eff_max-eff_min)/2.0;

    double statErr = 0.0;
    if (num > 0.0) statErr = eff * std::sqrt((1.0/num) + (1.0/denom));

    int n = out.pts->GetN();
    out.pts->SetPoint(n, x[b], eff);
    double errx = (b < xerrVec.size()) ? xerrVec[b] : 0.0;
    out.pts->SetPointError(n, errx, statErr);

    int m = out.band->GetN();
    out.band->SetPoint(m, x[b], eff);
    out.band->SetPointError(m, 0.0, thrErr);
  }
}

static void DrawGE1OnlyCanvas(
  TCanvas* c,
  const char* title,
  GraphPair& g,
  double yMax
){
  c->cd();
  g.band->SetTitle(title);
  g.band->GetXaxis()->SetTitle("Energy (MeV) (0.5 MeV bins)");
  g.band->GetYaxis()->SetTitle("Efficiency (%)");
  g.band->GetYaxis()->SetRangeUser(0, ResolveYMax(g, yMax));

  g.band->SetFillColorAlpha(kBlue, 0.20);
  g.band->SetFillStyle(1001);
  g.band->Draw("A3");

  g.pts->SetMarkerStyle(kFullDotLarge);
  g.pts->SetMarkerColor(kBlue);
  g.pts->SetLineColor(kBlue);
  g.pts->Draw("P");
}

static void DrawPerDetPads_GE1Only(
  TCanvas* c,
  const char* overallTitle,
  std::vector<GraphPair>& gDet,
  double yMax
){
  const int nDetectors = (int)gDet.size();
  if (nDetectors == 0) return;
  const int nColumns = (int)std::ceil(std::sqrt((double)nDetectors));
  const int nRows = (nDetectors + nColumns - 1) / nColumns;
  c->Divide(nColumns, nRows);

  for (int d=0; d<nDetectors; d++) {
    c->cd(d+1);

    gDet[d].band->SetTitle(Form("Liquid Scintillator %d; %s", d, overallTitle));
    gDet[d].band->GetXaxis()->SetTitle("Energy (MeV)");
    gDet[d].band->GetYaxis()->SetTitle("Efficiency (%)");
    gDet[d].band->GetYaxis()->SetRangeUser(0, ResolveYMax(gDet[d], yMax));

    gDet[d].band->SetFillColorAlpha(kBlue, 0.20);
    gDet[d].band->SetFillStyle(1001);
    gDet[d].band->Draw("A3");

    gDet[d].pts->SetMarkerStyle(kFullDotLarge);
    gDet[d].pts->SetMarkerColor(kBlue);
    gDet[d].pts->SetLineColor(kBlue);
    gDet[d].pts->Draw("P");
  }
}

void efficiency_total_perEne_1n(
  const char* inputFile = "bld/SECAR_sim_output.root"
)
{
  // -----------------------
  // USER SETTINGS
  // -----------------------
  const bool requireExactlyOneNeutron = true;  // safest for a true 1n dataset
  // Set either value > 0 to force that y-axis maximum (in percent).
  // Set it to 0 or a negative value to autoscale from the plotted data.
  // In auto mode, every per-detector graph is scaled independently.
  const double yMaxTotal = 0.0;
  const double yMaxDet   = 0.0;

  // One threshold value per detector. Add/remove entries as needed; all three
  // vectors must have the same size.
  const double heavy_thr  = 2.0;
  const std::vector<double> thr_min = {
    0.35, 0.35, 0.35, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5
  };
  const std::vector<double> thr_max = {
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
  };
  const std::vector<double> thr_nom = {
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8
  };
  if (thr_nom.empty() || thr_min.size() != thr_nom.size() ||
      thr_max.size() != thr_nom.size()) {
    std::cerr << "ERROR: threshold vectors must be non-empty and equal-sized.\n";
    return;
  }
  const int nDetectors = (int)thr_nom.size();

  // Energy binning
  // piecewise: 0.5 MeV bins up to 10 MeV, 2.0 MeV bins above 10 MeV
  const double Emin0   = 0.0;
  const double splitE  = 9.0;
  const double stepLow = 0.5;
  const double stepHigh= 2.0;

  // We'll build bins after we know the range of neutron energies for events kept
  std::vector<double> energyVec;                 // bin centers
  std::vector<double> xerrVec;                   // per-bin half-widths
  std::vector<double> binEdges;                  // edges: size = nbins+1

  auto buildBins = [&](double Emax)->void {
    binEdges.clear();
    energyVec.clear();
    xerrVec.clear();

    double edge = Emin0;
    binEdges.push_back(edge);
    // build bins until we cover up to Emax (include last bin containing Emax)
    while (edge < Emax) {
      double step = (edge >= splitE) ? stepHigh : stepLow;
      double next = edge + step;
      binEdges.push_back(next);
      energyVec.push_back(0.5*(edge+next));
      xerrVec.push_back((next-edge)/2.0);
      edge = next;
    }
  };

  auto getBin = [&](double E)->int {
    if (E < Emin0) return -1;
    if (binEdges.empty()) return -1;
    auto it = std::upper_bound(binEdges.begin(), binEdges.end(), E);
    if (it == binEdges.begin()) return -1;
    if (it == binEdges.end()) {
      // E is >= last edge: assign to last bin
      if (binEdges.size() < 2) return -1;
      return (int)binEdges.size() - 2;
    }
    int idx = (int)(it - binEdges.begin()) - 1;
    if (idx < 0) return -1;
    return idx;
  };

  // File / trees
  TFile* file = TFile::Open(inputFile, "READ");
  if (!file || file->IsZombie()) { std::cerr << "ERROR: Could not open ROOT file.\n"; return; }

  TTree* tTar = dynamic_cast<TTree*>(file->Get("TargetTrack"));
  TTree* tLS  = dynamic_cast<TTree*>(file->Get("LSHitTrack"));
  if (!tTar) { std::cerr << "ERROR: Missing TargetTrack\n"; file->Close(); return; }
  if (!tLS)  { std::cerr << "ERROR: Missing LSHitTrack\n";  file->Close(); return; }

  // Bind TargetTrack
  Double_t tar_eventID_d=0, tar_trackID_d=0, tar_isNeutron_d=0, tar_ekin_d=0;
  tTar->SetBranchAddress("eventID",   &tar_eventID_d);
  tTar->SetBranchAddress("trackID",   &tar_trackID_d);
  tTar->SetBranchAddress("isNeutron", &tar_isNeutron_d);
  tTar->SetBranchAddress("ekin_tar",  &tar_ekin_d);

  // Bind LSHitTrack
  Double_t ls_eventID_d=0, ls_detID_d=0, ls_parentID_d=0;
  Double_t ls_edep_d=0, ls_Zdep_d=0, ls_Adep_d=0;

  tLS->SetBranchAddress("eventID",   &ls_eventID_d);
  tLS->SetBranchAddress("detID",     &ls_detID_d);
  tLS->SetBranchAddress("parentID",  &ls_parentID_d);
  tLS->SetBranchAddress("edep_sum",  &ls_edep_d);
  tLS->SetBranchAddress("Z_dep",     &ls_Zdep_d);
  tLS->SetBranchAddress("A_dep",     &ls_Adep_d);

  // ------------------------------------------------------------
  // 1) Collect exactly one neutron per event: (event -> trackID, E)
  // ------------------------------------------------------------
  std::unordered_map<int, std::vector<int>> eventToNeutTracks;
  std::unordered_map<long long, double> neutronE;
  neutronE.reserve((size_t)tTar->GetEntries());

  const Long64_t nTar = tTar->GetEntries();
  for (Long64_t i=0; i<nTar; i++) {
    tTar->GetEntry(i);
    if ((int)std::llround(tar_isNeutron_d) != 1) continue;

    const int eventID = (int)std::llround(tar_eventID_d);
    const int trackID = (int)std::llround(tar_trackID_d);

    eventToNeutTracks[eventID].push_back(trackID);
    neutronE[Key32(eventID, trackID)] = (double)tar_ekin_d;
  }

  // Enforce 1 neutron/event (optionally)
  std::unordered_map<int, int> eventToOneTrack;
  eventToOneTrack.reserve(eventToNeutTracks.size());

  long long nEventsKept = 0, nEventsDropped = 0;
  for (auto& kv : eventToNeutTracks) {
    auto& v = kv.second;
    std::sort(v.begin(), v.end());

    if (requireExactlyOneNeutron) {
      if ((int)v.size() != 1) { nEventsDropped++; continue; }
      eventToOneTrack[kv.first] = v[0];
      nEventsKept++;
    } else {
      // if there are multiple neutrons, just take the first (lowest trackID)
      if (v.empty()) { nEventsDropped++; continue; }
      eventToOneTrack[kv.first] = v[0];
      nEventsKept++;
    }
  }

  std::cout << "\nTarget neutrons summary:\n";
  std::cout << "  Events kept:    " << nEventsKept << "\n";
  std::cout << "  Events dropped: " << nEventsDropped << "\n\n";

  // ------------------------------------------------------------
  // 2) Determine detection per event for TOTAL and for each detector
  //    We mark event detected if ANY LS hit passes threshold for that neutron track.
  // ------------------------------------------------------------
  // We'll size these after building the variable-width bins (based on event energies)
  std::vector<int> denomEvents;

  std::vector<int> numTot_nom, numTot_max, numTot_min;
  std::vector<std::vector<int>> numDet_nom(nDetectors), numDet_max(nDetectors),
                                numDet_min(nDetectors);

  // For integrated efficiencies
  long long N_events_total = 0;
  long long N_detTot_nom = 0;
  std::vector<long long> N_detDet_nom(nDetectors, 0);

  // Build bins based on actual event neutron energies (use events kept)
  double maxE = Emin0;
  for (const auto& kv : eventToOneTrack) {
    const int eventID = kv.first;
    const int trackID = kv.second;
    auto itE = neutronE.find(Key32(eventID, trackID));
    if (itE == neutronE.end()) continue;
    if (itE->second > maxE) maxE = itE->second;
  }
  if (maxE <= Emin0) { std::cerr << "ERROR: No valid event energies to bin.\n"; file->Close(); return; }

  buildBins(maxE);
  int nBins = (int)energyVec.size();
  if (nBins <= 0) { std::cerr << "ERROR: No bins created.\n"; file->Close(); return; }

  // allocate per-bin counters
  denomEvents.assign(nBins, 0);
  numTot_nom.assign(nBins, 0);
  numTot_max.assign(nBins, 0);
  numTot_min.assign(nBins, 0);
  for (int d=0; d<nDetectors; d++){
    numDet_nom[d].assign(nBins,0);
    numDet_max[d].assign(nBins,0);
    numDet_min[d].assign(nBins,0);
  }

  // Build a quick lookup: event -> (trackID, bin)
  struct EvInfo { int trackID; int bin; };
  std::unordered_map<int, EvInfo> evInfo;
  evInfo.reserve(eventToOneTrack.size());

  for (const auto& kv : eventToOneTrack) {
    const int eventID = kv.first;
    const int trackID = kv.second;
    auto itE = neutronE.find(Key32(eventID, trackID));
    if (itE == neutronE.end()) continue;

    int b = getBin(itE->second);
    if (b < 0) continue;

    evInfo[eventID] = {trackID, b};
    denomEvents[b]++;
    N_events_total++;
  }

  // Flags per event for detection (nom/max/min), total and per det
  struct Flags {
    bool tot_nom=false, tot_max=false, tot_min=false;
    std::vector<bool> det_nom, det_max, det_min;
    explicit Flags(int n = 0)
      : det_nom(n, false), det_max(n, false), det_min(n, false) {}
  };
  std::unordered_map<int, Flags> flags;
  flags.reserve(evInfo.size());

  const Long64_t nLS = tLS->GetEntries();
  for (Long64_t i=0; i<nLS; i++) {
    tLS->GetEntry(i);

    const int eventID = (int)std::llround(ls_eventID_d);
    auto itEv = evInfo.find(eventID);
    if (itEv == evInfo.end()) continue;

    const int detID = (int)std::llround(ls_detID_d);
    if (detID < 0 || detID >= nDetectors) continue;

    const int parentID = (int)std::llround(ls_parentID_d);
    if (parentID != itEv->second.trackID) continue; // correlate hit to the one neutron

    const double edep = (double)ls_edep_d;
    const int Zdep    = (int)std::llround(ls_Zdep_d);
    const int Adep    = (int)std::llround(ls_Adep_d);

    const bool pass_nom = PassThreshold(edep, Zdep, Adep, thr_nom[detID], heavy_thr);
    const bool pass_max = PassThreshold(edep, Zdep, Adep, thr_max[detID], heavy_thr);
    const bool pass_min = PassThreshold(edep, Zdep, Adep, thr_min[detID], heavy_thr);

    auto inserted = flags.emplace(eventID, Flags(nDetectors));
    auto& f = inserted.first->second;
    if (pass_nom) { f.tot_nom=true; f.det_nom[detID]=true; }
    if (pass_max) { f.tot_max=true; f.det_max[detID]=true; }
    if (pass_min) { f.tot_min=true; f.det_min[detID]=true; }
  }

  // Convert flags -> binned numerators
  for (const auto& kv : evInfo) {
    const int eventID = kv.first;
    const int b = kv.second.bin;
    auto itF = flags.find(eventID);
    if (itF == flags.end()) continue;

    const auto& f = itF->second;

    if (f.tot_nom) { numTot_nom[b]++; N_detTot_nom++; }
    if (f.tot_max) numTot_max[b]++;
    if (f.tot_min) numTot_min[b]++;

    for(int d=0; d<nDetectors; d++){
      if (f.det_nom[d]) {numDet_nom[d][b]++; N_detDet_nom[d]++; }
      if (f.det_max[d]) numDet_max[d][b]++;
      if (f.det_min[d]) numDet_min[d][b]++;
    }
  }

  // ------------------------------------------------------------
  // Multiplicity summary (integrated over all bins) for nominal thresholds
  // Count how many detectors fired per event (should be mostly 1 for 1n datasets)
  long long N_mult1_nom = 0;
  long long N_mult_gt1_nom = 0;
  long long N_mult_gt2_nom = 0;
  for (const auto& kv : evInfo) {
    const int eventID = kv.first;
    auto itF = flags.find(eventID);
    if (itF == flags.end()) continue;
    const auto& f = itF->second;
    int cnt = 0;
    for (int d=0; d<nDetectors; d++) if (f.det_nom[d]) cnt++;
    if (cnt == 1) N_mult1_nom++;
    if (cnt > 1)  N_mult_gt1_nom++;
    if (cnt > 2)  N_mult_gt2_nom++;
  }

  // Multiplicity summary for max/min thresholds
  long long N_mult1_max = 0;
  long long N_mult_gt1_max = 0;
  long long N_mult_gt2_max = 0;
  long long N_mult1_min = 0;
  long long N_mult_gt1_min = 0;
  long long N_mult_gt2_min = 0;
  for (const auto& kv : evInfo) {
    const int eventID = kv.first;
    auto itF = flags.find(eventID);
    if (itF == flags.end()) continue;
    const auto& f = itF->second;
    int cnt_max = 0, cnt_min = 0;
    for (int d=0; d<nDetectors; d++) {
      if (f.det_max[d]) cnt_max++;
      if (f.det_min[d]) cnt_min++;
    }
    if (cnt_max == 1) N_mult1_max++;
    if (cnt_max > 1)  N_mult_gt1_max++;
    if (cnt_max > 2)  N_mult_gt2_max++;
    if (cnt_min == 1) N_mult1_min++;
    if (cnt_min > 1)  N_mult_gt1_min++;
    if (cnt_min > 2)  N_mult_gt2_min++;
  }

  // ------------------------------------------------------------
  // 3) Build graphs + plot
  // ------------------------------------------------------------
  GraphPair gTot = MakePair();
  FillEventGraph(energyVec, denomEvents, numTot_nom, numTot_max, numTot_min, xerrVec, gTot);

  std::vector<GraphPair> gDet(nDetectors);
  for(int d=0; d<nDetectors; d++){
    gDet[d] = MakePair();
    FillEventGraph(energyVec, denomEvents, numDet_nom[d], numDet_max[d], numDet_min[d], xerrVec, gDet[d]);
  }

  TCanvas* cTot = new TCanvas("cTot_1n", "TOTAL (any LS) - 1n/event", 1100, 800);
  DrawGE1OnlyCanvas(cTot, "Cumulative efficiency", gTot, yMaxTotal);

  TCanvas* cDet = new TCanvas("cDet_1n", "Per-detector - 1n/event", 1200, 800);
  DrawPerDetPads_GE1Only(cDet, "1 neutron/event", gDet, yMaxDet);

  // ------------------------------------------------------------
  // 4) Print integrated efficiencies (nominal, max, min)
  // ------------------------------------------------------------
  auto sumVec = [&](const std::vector<int>& v)->long long { long long s=0; for(int x: v) s += x; return s; };

  // totals already tracked for nominal in N_detTot_nom / N_detDet_nom
  const long long N_detTot_max = sumVec(numTot_max);
  std::vector<long long> N_detDet_max(nDetectors, 0);
  for (int d=0; d<nDetectors; d++) N_detDet_max[d] = sumVec(numDet_max[d]);

  const long long N_detTot_min = sumVec(numTot_min);
  std::vector<long long> N_detDet_min(nDetectors, 0);
  for (int d=0; d<nDetectors; d++) N_detDet_min[d] = sumVec(numDet_min[d]);

  auto safeEff = [&](long long n)->double {
    if (N_events_total<=0) return 0.0;
    return 100.0 * (double)n / (double)N_events_total;
  };

  std::cout << "\n==================== FINAL EVENT EFFICIENCIES (1 neutron/event) ====================\n";
  std::cout << "Total events used: " << N_events_total << "\n\n";

  std::cout << "NOMINAL thresholds:\n";
  std::cout << "  TOTAL (any LS) : " << N_detTot_nom << "  =>  " << safeEff(N_detTot_nom) << " %\n";
  // for(int d=0; d<4; d++){
  //   std::cout << "  LS" << d << "           : " << N_detDet_nom[d] << "  =>  " << safeEff(N_detDet_nom[d]) << " %\n";
  // }
  std::cout << "\nMultiplicity summary (nominal thresholds):\n";
  std::cout << "  Mult = 1 : " << N_mult1_nom << "  =>  " << safeEff(N_mult1_nom) << " %\n";
  std::cout << "  Mult > 1 : " << N_mult_gt1_nom << "  =>  " << safeEff(N_mult_gt1_nom) << " %\n";
  std::cout << "  Mult > 2 : " << N_mult_gt2_nom << "  =>  " << safeEff(N_mult_gt2_nom) << " %\n\n";

  std::cout << "MAX thresholds:\n";
  std::cout << "  TOTAL (any LS) : " << N_detTot_max << "  =>  " << safeEff(N_detTot_max) << " %\n";
  // for(int d=0; d<4; d++){
  //   std::cout << "  LS" << d << "           : " << N_detDet_max[d] << "  =>  " << safeEff(N_detDet_max[d]) << " %\n";
  // }
  std::cout << "\nMultiplicity summary (MAX thresholds):\n";
  std::cout << "  Mult = 1 : " << N_mult1_max << "  =>  " << safeEff(N_mult1_max) << " %\n";
  std::cout << "  Mult > 1 : " << N_mult_gt1_max << "  =>  " << safeEff(N_mult_gt1_max) << " %\n";
  std::cout << "  Mult > 2 : " << N_mult_gt2_max << "  =>  " << safeEff(N_mult_gt2_max) << " %\n";
  std::cout << "\nMIN thresholds:\n";
  std::cout << "  TOTAL (any LS) : " << N_detTot_min << "  =>  " << safeEff(N_detTot_min) << " %\n";
  // for(int d=0; d<4; d++){
  //   std::cout << "  LS" << d << "           : " << N_detDet_min[d] << "  =>  " << safeEff(N_detDet_min[d]) << " %\n";
  // }
  std::cout << "\nMultiplicity summary (MIN thresholds):\n";
  std::cout << "  Mult = 1 : " << N_mult1_min << "  =>  " << safeEff(N_mult1_min) << " %\n";
  std::cout << "  Mult > 1 : " << N_mult_gt1_min << "  =>  " << safeEff(N_mult_gt1_min) << " %\n";
  std::cout << "  Mult > 2 : " << N_mult_gt2_min << "  =>  " << safeEff(N_mult_gt2_min) << " %\n";

  std::cout << "=====================================================================================" << "\n\n";

  file->Close();
}
