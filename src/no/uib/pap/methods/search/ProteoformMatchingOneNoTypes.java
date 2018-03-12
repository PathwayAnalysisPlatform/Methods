package no.uib.pap.methods.search;

import no.uib.pap.model.Proteoform;
import org.apache.commons.lang3.tuple.Pair;

public class ProteoformMatchingOneNoTypes extends ProteoformMatching {
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

        if (rP.getPtms().size() == 0) {
            return true;
        }

        // At least one of the reference ptms should be in the input
        for (Pair<String, Long> rPtm : rP.getPtms()) {
            for (Pair<String, Long> iPtm : iP.getPtms()) {
                // If corrdinate matches
                if (matches(rPtm.getValue(), iPtm.getValue(), margin)) {
                    return true;
                }
            }
        }

        return false;
    }
}