package no.uib.pap.methods.analysis.ora;

import com.google.common.collect.ImmutableMap;
import no.uib.pap.methods.search.Search;
import no.uib.pap.methods.search.SearchResult;
import no.uib.pap.model.InputType;
import no.uib.pap.model.MessageStatus;
import no.uib.pap.model.Pathway;
import org.apache.commons.math3.distribution.BinomialDistribution;

import java.util.Comparator;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;

public class Analysis {

    /**
     * Performs over representation analysis on the hit pathways by the search.
     *
     * @param searchResult   Pathway instances with the found entities and the counts. Results of the analysis go inside this instances.
     * @param populationSize Total number of proteins(counting isoform) or proteoforms in Reactome
     * @return Error or success status messages
     */
    public static AnalysisResult analysis(SearchResult searchResult, int populationSize) {

        // Traverse all the pathways
        int percentage = 0;
        int processed = 0;
        for (Pathway pathway : searchResult.getHitPathways()) {

            // Calculate proteoformSet and iReactions ratio
            pathway.setEntitiesRatio(
                    (double) pathway.getEntitiesFound().size() / (double) pathway.getNumEntitiesTotal());
            pathway.setReactionsRatio(
                    (double) pathway.getReactionsFound().size() / (double) pathway.getNumReactionsTotal());

            // Calculate the proteoformSet pvalue
            int k = pathway.getEntitiesFound().size(); // Sucessful trials: Entities found participating in the pathway
            double p = pathway.getNumEntitiesTotal() / (double) populationSize; // Probability of sucess in each trial: Entities in the pathway / All possible entities

            BinomialDistribution binomialDistribution = new BinomialDistribution(searchResult.getHitProteins().size(), p); // Given n trials with probability p of success
            pathway.setpValue(1 - binomialDistribution.cumulativeProbability(k - 1)); // Probability of k or more successful trials

            processed++;
            int newPercentage = processed * 100 / searchResult.getHitPathways().size();
            if (newPercentage > percentage + 2) {
                System.out.print(newPercentage + "% ");
                percentage = newPercentage;
            }
        }
        System.out.println("\n");

        adjustPValues(searchResult.getHitPathways());

        return new AnalysisResult(searchResult.getHitPathways(), new MessageStatus("Sucess", 0, 0, "", ""));
    }

    public static int getPopulationSize(InputType inputType, int totalProteins, int totalProteoforms) {
        switch (inputType) {
            case GENE:
            case GENES:
            case ENSEMBL:
            case ENSEMBLS:
            case UNIPROT:
            case UNIPROTS:
            case RSID:
            case RSIDS:
            case CHRBP:
            case CHRBPS:
            case VCF:
            case PEPTIDE:
            case PEPTIDES:
                return totalProteins;
            case PROTEOFORM:
            case PROTEOFORMS:
            case MODIFIEDPEPTIDE:
            case MODIFIEDPEPTIDES:
                return totalProteoforms;
            default:
                return 0;
        }
    }

    /**
     * Benjamini-Hochberg adjustment for FDR at 0.05%
     */
    private static void adjustPValues(Set<Pathway> hitPathways) {

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

        for (Pathway pathway : hitPathways) {
            sortedPathways.add(pathway);
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
