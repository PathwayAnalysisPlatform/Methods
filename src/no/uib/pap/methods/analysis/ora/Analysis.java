package no.uib.pap.methods.analysis.ora;

import com.google.common.collect.ImmutableMap;
import no.uib.pap.model.MessageStatus;
import no.uib.pap.model.Pathway;
import no.uib.pap.model.ProteoformFormat;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.BinomialDistribution;

import java.text.ParseException;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.TreeSet;

public class Analysis {

    public static Pair<TreeSet<Pathway>, MessageStatus> analysis(
            ImmutableMap<String, Pathway> iPathways,
            int populationSize,
            List<String[]> searchResult) {

        MessageStatus status = null;

        HashSet<String> hitProteins = new HashSet<String>();
        HashSet<String> hitPathways = new HashSet<String>();


        // Get the set of hit pathways and calculate the found entities for each pathway
        for (String[] record : searchResult) {
            String entity = record[1];
            String reactionStId = record[3];
            String pathwayStId = record[5];

            hitProteins.add(entity);
            hitPathways.add(pathwayStId);

            // Add current protein to the found proteoformSet of the pathway
            Pathway pathway = iPathways.get(pathwayStId);
            pathway.getReactionsFound().add(reactionStId);
            try {
                pathway.getEntitiesFound().add(ProteoformFormat.SIMPLE.getProteoform(entity, 0));
            } catch (ParseException e) {
                e.printStackTrace(); // TODO Send correct error
            }
        }

        // Traverse all the iPathways
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

        }

        return new MutablePair<TreeSet<Pathway>, MessageStatus>(adjustPValues(iPathways, hitPathways), new MessageStatus("Sucess", 0, 0, "", ""));
    }

    /**
     * Benjamini-Hochberg adjustment for FDR at 0.05%
     */
    private static TreeSet<Pathway> adjustPValues(ImmutableMap<String, Pathway> iPathways, HashSet<String> hitPathways) {

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

        return sortedPathways;
    }
}
