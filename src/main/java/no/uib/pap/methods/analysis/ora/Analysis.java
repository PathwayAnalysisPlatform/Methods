package no.uib.pap.methods.analysis.ora;

import com.google.common.collect.ImmutableMap;
import no.uib.pap.model.MessageStatus;
import no.uib.pap.model.Pathway;
import org.apache.commons.math3.distribution.BinomialDistribution;

import java.util.Comparator;
import java.util.HashSet;
import java.util.TreeSet;

public class Analysis {

    /**
     * Performs over representation analysis on the hit pathways by the search.
     * @param iPathways Pathway instances with the found entities and the counts. Results of the analysis go inside this instances.
     * @param populationSize Total number of proteins(counting isoform) or proteoforms in Reactome
     * @param hitProteins Participant proteins in the hit pathways
     * @param hitPathways Selected/hit pathways to be analysed
     * @return Error or success status messages
     */
    public static MessageStatus analysis(
            ImmutableMap<String, Pathway> iPathways,
            int populationSize,
            TreeSet<String> hitProteins,
            HashSet<String> hitPathways) {

        // Traverse all the iPathways
        int percentage = 0;
        int processed = 0;
        for (String stId : hitPathways) {

            Pathway pathway = iPathways.get(stId);

            // Calculate proteoformSet and iReactions ratio
            pathway.setEntitiesRatio(
                    (double) pathway.getEntitiesFound().size() / (double) pathway.getNumEntitiesTotal());
            pathway.setReactionsRatio(
                    (double) pathway.getReactionsFound().size() / (double) pathway.getNumReactionsTotal());

            // Calculate the proteoformSet pvalue
            int k = pathway.getEntitiesFound().size(); // Sucessful trials: Entities found participating in the pathway
            double p = pathway.getNumEntitiesTotal() / (double) populationSize; // Probability of sucess in each trial: The entity is a participant in the pathway

            BinomialDistribution binomialDistribution = new BinomialDistribution(hitProteins.size(), p); // Given n trials with probability p of success
            pathway.setpValue(binomialDistribution.probability(k)); // Probability of k successful trials

            processed++;
            int newPercentage = processed * 100 / hitPathways.size();
            if(newPercentage > percentage + 2){
                System.out.print(newPercentage + "% ");
                percentage = newPercentage;
            }
        }
        System.out.println("\n");

        adjustPValues(iPathways, hitPathways);

        return new MessageStatus("Sucess", 0, 0, "", "");
    }

    /**
     * Benjamini-Hochberg adjustment for FDR at 0.05%
     */
    private static void adjustPValues(ImmutableMap<String, Pathway> iPathways, HashSet<String> hitPathways) {

        // Sort iPathways by pValue
        Comparator<Pathway> comparator = new Comparator<Pathway>() {
            public int compare(Pathway x, Pathway y) {

                if (x.equals(y))
                    return 0;

                if (x.getPValue() != y.getPValue()) {
                    return Double.compare(x.getPValue(), y.getPValue());
                }

                // First by displayName
                if (!x.getDisplayName().equals(y.getDisplayName())) {
                    return x.getDisplayName().compareTo(y.getDisplayName());
                }

                // Second by stId
                if (!x.getStId().equals(y.getStId())) {
                    return x.getStId().compareTo(y.getStId());
                }

                return 0;
            }
        };

        TreeSet<Pathway> sortedPathways = new TreeSet<Pathway>(comparator);

        for (String stId : hitPathways) {
            sortedPathways.add(iPathways.get(stId));
        }
        // System.out.println("The number of pathways to be analysed is: " +
        // sortedPathways.size());
        // Count number of iPathways with p-Values less than 0.05
        double n = 0;
        for (Pathway pathway : sortedPathways) {
            if (pathway.getPValue() < 0.05) {
                n++;
            } else {
                break;
            }
        }

        double rank = 1;
        for (Pathway pathway : sortedPathways) {
            double newPValue = pathway.getPValue() * n;
            newPValue /= rank;
            pathway.setEntitiesFDR(newPValue);
            rank++;
        }
        System.out.println("The number of analysed pathways is: " + sortedPathways.size());
    }
}
