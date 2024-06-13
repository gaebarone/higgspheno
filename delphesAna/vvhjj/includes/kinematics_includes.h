#ifndef KINEMATICS_INCLUDE_H
#define KINEMATICS_INCLUDE_H



double deltaPhi(TLorentzVector &lorentzvector1, TLorentzVector &lorentzvector2) {

    // return (lorentzvector1.Phi() > lorentzvector2.Phi() ? -1:+1)*TMath::Abs(lorentzvector2.Phi() - lorentzvector1.Phi());
    return (lorentzvector1.Phi() - lorentzvector2.Phi());
}

double deltaEta(TLorentzVector &lorentzvector1, TLorentzVector &lorentzvector2) {

    return (lorentzvector1.Eta() - lorentzvector2.Eta());

}

double deltaR(TLorentzVector &lorentzvector1, TLorentzVector &lorentzvector2) {

    double deltaPhi = lorentzvector1.Phi() - lorentzvector2.Phi();
    double deltaEta = lorentzvector1.Eta() - lorentzvector2.Eta();
    return std::sqrt(deltaPhi**2 + deltaEta**2);
    
}

double massTransverse(TLorentzVector &lorentzvector1, TLorentzVector &lorentzvector2) {


    double pT_lepton = lorentzvector1.Pt()
    double pT_miss = lorentzvector2.Pt()
    double deltaPhi = deltaPhi(lorentzvector1, lorentzvector2)
    return std::sqrt(2 * pT_lepton * pT_miss * (1 - std::cos(deltaPhi)));

}

