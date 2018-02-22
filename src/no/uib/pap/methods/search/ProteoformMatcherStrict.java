package no.uib.pap.methods.search;

import no.uib.pap.model.Proteoform;
import org.apache.commons.lang3.tuple.MutablePair;

import java.util.Map;

public class ProteoformMatcherStrict extends ProteoformMatcher {

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

        if (rP.getPtms().size() != iP.getPtms().size()) {
            return false;
        }

        // All the reference PTMs should be exactly in the input
        for (Map.Entry<String, Long> rPtm : rP.getPtms()) {
            if (!iP.getPtms().contains(new MutablePair<>(rPtm.getKey(), rPtm.getValue()))) {
                return false;
            }
        }

        // All the input PTMs should be exactly in the reference
        for (Map.Entry<String, Long> iPtm : iP.getPtms()) {
            if (!rP.getPtms().contains(new MutablePair<>(iPtm.getKey(), iPtm.getValue()))) {
                return false;
            }
        }

        return true;
    }

}
