package no.uib.pap.methods.search;

import no.uib.pap.model.Proteoform;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;

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
        for (Pair<String, Long> rPtm : rP.getPtms()) {
            if (!iP.getPtms().contains(new MutablePair<>(rPtm.getKey(), rPtm.getValue()))) {
                boolean anyMatches = false;
                for (Pair<String, Long> iPtm : iP.getPtms()) {
                    if (rPtm.getLeft().equals(iPtm.getLeft())) {
                        if (matches(rPtm.getRight(), iPtm.getRight(), margin)) {
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
