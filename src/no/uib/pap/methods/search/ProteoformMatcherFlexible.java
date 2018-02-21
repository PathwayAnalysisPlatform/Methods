package no.uib.pap.methods.search;

import no.uib.pap.model.Proteoform;

import java.util.Map;

public class ProteoformMatcherFlexible extends ProteoformMatcher {

	public boolean matches(Long iC, Long rC, Long margin){
        if(iC != null){ if(iC == -1L) iC = null; }
        if(rC != null){ if(rC == -1L) rC = null; }
        if(iC != null && rC != null){
            if(iC != rC){
                if(Math.abs(iC-rC) > margin){
                    return false;
                }
            }
        }
        return true;
    }
	
    @Override
    public Boolean matches(Proteoform iP, Proteoform rP, Long margin) {

        // Check the uniprot accession, including the isoform matches
        if (iP.getUniProtAcc() == null) {
            throw new IllegalArgumentException();
        }

        if (rP.getUniProtAcc() == null) {
            throw new IllegalArgumentException();
        }

        if (!iP.getUniProtAcc().equals(rP.getUniProtAcc())) {
            return false;
        }

        if (!matches(iP.getStartCoordinate(), rP.getStartCoordinate(), margin)) {
            return false;
        }

        if (!matches(iP.getEndCoordinate(), rP.getEndCoordinate(), margin)) {
            return false;
        }

        // All the reference PTMs should be in the input
        for (Map.Entry<String, Long> rPtm : rP.getPtms().entries()) {
            if (!iP.getPtms().containsEntry(rPtm.getKey(), rPtm.getValue())) {
                boolean anyMatches = false;
                for (Map.Entry<String, Long> iPtm : iP.getPtms().entries()) {
                    if (rPtm.getKey().equals(iPtm.getKey())) {
                        if (matches(rPtm.getValue(), iPtm.getValue(), margin)) {
                            anyMatches = true;
                            break;
                        }
                    }
                }
                if (!anyMatches) {
                    return false;
                }
            }
        }

        return true;
    }

}
