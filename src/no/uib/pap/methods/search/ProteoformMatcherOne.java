package no.uib.pap.methods.search;

import no.uib.pap.model.Proteoform;

import java.util.Map;

public class ProteoformMatcherOne extends ProteoformMatcher {
    @Override
    public Boolean matches(Proteoform iP, Proteoform rP, Long margin) {

        // Check the uniprot accession, including the isoform matches
        if(iP.getUniProtAcc() == null){
            throw new IllegalArgumentException();
        }

        if(rP.getUniProtAcc() == null){
            throw new IllegalArgumentException();
        }

        if(!iP.getUniProtAcc().equals(rP.getUniProtAcc())){
            return false;
        }

        if(!matches(iP.getStartCoordinate(), rP.getStartCoordinate(), margin)){
            return false;
        }

        if(!matches(iP.getEndCoordinate(), rP.getEndCoordinate(), margin)){
            return false;
        }

        if(rP.getPtms().entries().size() == 0){
            return true;
        }

        // At least one of the reference ptms should be in the input
        for(Map.Entry<String, Long> rPtm : rP.getPtms().entries()){
            if(iP.getPtms().containsEntry(rPtm.getKey(), rPtm.getValue())){
                return true;
            }
            for (Map.Entry<String, Long> iPtm : iP.getPtms().entries()) {
                if (rPtm.getKey().equals(iPtm.getKey())) {
                    if (matches(rPtm.getValue(), iPtm.getValue(), margin)) {
                        return true;
                    }
                }
            }
        }

        return false;
    }
}
