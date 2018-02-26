package no.uib.pap.methods.search;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSetMultimap;
import no.uib.pap.model.*;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;

import java.lang.Error;
import java.text.ParseException;
import java.util.*;

import static no.uib.pap.methods.search.PeptideMatcher.*;
import static no.uib.pap.model.Error.ERROR_INITIALIZING_PEPTIDE_MAPPER;
import static no.uib.pap.model.InputPatterns.*;
import static no.uib.pap.model.Warning.*;

/**
 * Methods to get the reactions and pathways using a list of entities of the accepted input types.
 * <p>
 * <p>These methods must fill the number of reactions and entities found in each pathway in the structure passed as parameter.
 * They also must fill in the set for the hit proteins and hit pathways.</p>
 */
public class Search {

    public static HashSet<String> hitProteins = new HashSet<>(); // These are in the reference data
    static final String fasta = "uniprot-all.fasta";

    // Fills the hitProteins set to call the next method
    public static Pair<List<String[]>, MessageStatus> searchWithUniProt(
            Collection<String> input,
            ImmutableMap<String, String> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, String> imapProteinsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            HashSet<String> hitProteins,
            HashSet<String> hitPathways) {

       List<String[]> result = new ArrayList<String[]>();

        for (String protein : input) {
            protein = protein.trim();

            for (String reaction : imapProteinsToReactions.get(protein)) {
                hitProteins.add(protein);

                for (String pathwayStId : imapReactionsToPathways.get(reaction)) {
                    hitPathways.add(pathwayStId);
                    Pathway pathway = iPathways.get(pathwayStId);
                    pathway.getReactionsFound().add(reaction);
                    try {
                        pathway.getEntitiesFound().add(ProteoformFormat.SIMPLE.getProteoform(protein, 0));
                    } catch (ParseException e) {
                        return new MutablePair<List<String[]>, MessageStatus>(
                                result,
                                new MessageStatus(
                                        "Failed",
                                        no.uib.pap.model.Error.INPUT_PARSING_ERROR.getCode(),
                                        no.uib.pap.model.Error.INPUT_PARSING_ERROR.getCode(),
                                        no.uib.pap.model.Error.INPUT_PARSING_ERROR.getMessage(),
                                        ""));
                    }
                    if (topLevelPathways && imapPathwaysToTopLevelPathways.get(pathwayStId).size() > 0) {
                        for (String topLevelPathway : imapPathwaysToTopLevelPathways.get(pathwayStId)) {
                            String[] values = {
                                    "",
                                    protein,
                                    reaction,
                                    iReactions.get(reaction),
                                    pathwayStId,
                                    iPathways.get(pathwayStId).getDisplayName(),
                                    topLevelPathway,
                                    iPathways.get(topLevelPathway).getDisplayName()
                            };
                            result.add(values);
                        }
                    } else {
                        String[] values = {
                                "",
                                protein,
                                reaction,
                                iReactions.get(reaction),
                                pathwayStId,
                                iPathways.get(pathwayStId).getDisplayName()
                        };
                        result.add(values);
                    }
                }
            }
        }

        System.out.println("Requested " + hitProteins.size() + " proteins.");

        MessageStatus status = null;
        status = new MessageStatus("Sucess", 0, 0, "", "");
        return new MutablePair<>(result, status);
    }


    public static Pair<List<String[]>, MessageStatus> searchWithGene(
            List<String> input,
            ImmutableMap<String, String> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, String> imapGenesToProteins,
            ImmutableSetMultimap<String, String> imapProteinsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            HashSet<String> hitProteins,
            HashSet<String> hitPathways) {

        List<String[]> result = new ArrayList<String[]>();

        for (String gene : input) {
            for (String protein : imapGenesToProteins.get(gene.trim())) {
                hitProteins.add(protein);
                for (String reaction : imapProteinsToReactions.get(protein)) {
                    for (String pathwayStId : imapReactionsToPathways.get(reaction)) {
                        hitPathways.add(pathwayStId);
                        Pathway pathway = iPathways.get(pathwayStId);
                        pathway.getReactionsFound().add(reaction);
                        try {
                            pathway.getEntitiesFound().add(ProteoformFormat.SIMPLE.getProteoform(protein, 0));
                        } catch (ParseException e) {
                            return new MutablePair<List<String[]>, MessageStatus>(
                                    result,
                                    new MessageStatus(
                                            "Failed",
                                            no.uib.pap.model.Error.INPUT_PARSING_ERROR.getCode(),
                                            no.uib.pap.model.Error.INPUT_PARSING_ERROR.getCode(),
                                            no.uib.pap.model.Error.INPUT_PARSING_ERROR.getMessage(),
                                            ""));
                        }
                        if (topLevelPathways && imapPathwaysToTopLevelPathways.get(pathwayStId).size() > 0) {
                            for (String topLevelPathway : imapPathwaysToTopLevelPathways.get(pathwayStId)) {
                                String[] values = {
                                        gene,
                                        protein,
                                        reaction,
                                        iReactions.get(reaction),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName(),
                                        topLevelPathway,
                                        iPathways.get(topLevelPathway).getDisplayName()
                                };
                                result.add(values);
                            }
                        } else {
                            String[] values = {
                                    gene,
                                    protein,
                                    reaction,
                                    iReactions.get(reaction),
                                    pathwayStId,
                                    iPathways.get(pathwayStId).getDisplayName()
                            };
                            result.add(values);
                        }
                    }
                }
            }
        }

        System.out.println("Requested " + hitProteins.size() + " proteins.");

        MessageStatus status = null;
        status = new MessageStatus("Sucess", 0, 0, "", "");
        return new MutablePair<>(result, status);
    }

    public static Pair<List<String[]>, MessageStatus> searchWithEnsembl(
            List<String> input,
            ImmutableMap<String, String> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, String> imapEnsemblToProteins,
            ImmutableSetMultimap<String, String> imapProteinsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            HashSet<String> hitProteins,
            HashSet<String> hitPathways) {

        List<String[]> result = new ArrayList<String[]>();

        for (String ensembl : input) {
            for (String protein : imapEnsemblToProteins.get(ensembl)) {
                for (String reaction : imapProteinsToReactions.get(protein)) {
                    for (String pathwayStId : imapReactionsToPathways.get(reaction)) {

                        if (topLevelPathways && imapPathwaysToTopLevelPathways.get(pathwayStId).size() > 0) {
                            for (String topLevelPathway : imapPathwaysToTopLevelPathways.get(pathwayStId)) {
                                String[] values = {
                                        ensembl,
                                        protein,
                                        reaction,
                                        iReactions.get(reaction),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName(),
                                        topLevelPathway,
                                        iPathways.get(topLevelPathway).getDisplayName()
                                };
                                result.add(values);
                            }
                        } else {
                            String[] values = {
                                    ensembl,
                                    protein,
                                    reaction,
                                    iReactions.get(reaction),
                                    pathwayStId,
                                    iPathways.get(pathwayStId).getDisplayName()
                            };
                            result.add(values);
                        }
                    }
                }
            }
        }

        System.out.println("Requested " + hitProteins.size() + " proteins.");

        MessageStatus status = null;
        status = new MessageStatus("Sucess", 0, 0, "", "");
        return new MutablePair<>(result, status);
    }

    public static Pair<List<String[]>, MessageStatus> searchWithRsId(
            List<String> input,
            ImmutableMap<String, String> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, String> imapRsIdsToProteins,
            ImmutableSetMultimap<String, String> imapProteinsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            HashSet<String> hitProteins,
            HashSet<String> hitPathways) {

        List<String[]> result = new ArrayList<String[]>();

        int row = 0;
        for (String rsid : input) {
            row++;
            if (rsid.isEmpty()) {
                sendWarning(EMPTY_ROW, row);
                continue;
            }
            if (!matches_Rsid(rsid)) {
                sendWarning(INVALID_ROW, row);
                continue;
            }
            for (String protein : imapRsIdsToProteins.get(rsid)) {
                for (String reaction : imapProteinsToReactions.get(protein)) {
                    for (String pathwayStId : imapReactionsToPathways.get(reaction)) {

                        if (topLevelPathways && imapPathwaysToTopLevelPathways.get(pathwayStId).size() > 0) {
                            for (String topLevelPathway : imapPathwaysToTopLevelPathways.get(pathwayStId)) {
                                String[] values = {
                                        rsid,
                                        protein,
                                        reaction,
                                        iReactions.get(reaction),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName(),
                                        topLevelPathway,
                                        iPathways.get(topLevelPathway).getDisplayName()
                                };
                                result.add(values);
                            }
                        } else {
                            String[] values = {
                                    rsid,
                                    protein,
                                    reaction,
                                    iReactions.get(reaction),
                                    pathwayStId,
                                    iPathways.get(pathwayStId).getDisplayName()
                            };
                            result.add(values);
                        }
                    }
                }
            }
        }

        System.out.println("Requested " + hitProteins.size() + " proteins.");

        MessageStatus status = null;
        status = new MessageStatus("Sucess", 0, 0, "", "");
        return new MutablePair<>(result, status);
    }

    public static Pair<List<String[]>, MessageStatus> searchWithChrBp(
            List<String> input,
            ImmutableMap<String, String> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, String> imapChrBpToProteins,
            ImmutableSetMultimap<String, String> imapProteinsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            HashSet<String> hitProteins,
            HashSet<String> hitPathways) {

        List<String[]> result = new ArrayList<String[]>();

        int row = 0;
        for (String line : input) {
            row++;
            if (line.isEmpty()) {
                sendWarning(EMPTY_ROW, row);
                continue;
            }
            if (!matches_ChrBp(line)) {
                sendWarning(INVALID_ROW, row);
                continue;
            }

            Snp snp = getSnpFromChrBp(line);
            for (String protein : imapChrBpToProteins.get(snp.getChr() + "_" + snp.getBp())) {
                for (String reaction : imapProteinsToReactions.get(protein)) {
                    for (String pathwayStId : imapReactionsToPathways.get(reaction)) {

                        if (topLevelPathways && imapPathwaysToTopLevelPathways.get(pathwayStId).size() > 0) {
                            for (String topLevelPathway : imapPathwaysToTopLevelPathways.get(pathwayStId)) {
                                String[] values = {
                                        snp.getChr() + "+" + snp.getBp(),
                                        protein,
                                        reaction,
                                        iReactions.get(reaction),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName(),
                                        topLevelPathway,
                                        iPathways.get(topLevelPathway).getDisplayName()
                                };
                                result.add(values);
                            }
                        } else {
                            String[] values = {
                                    snp.getChr() + "+" + snp.getBp(),
                                    protein,
                                    reaction,
                                    iReactions.get(reaction),
                                    pathwayStId,
                                    iPathways.get(pathwayStId).getDisplayName()
                            };
                            result.add(values);
                        }
                    }
                }
            }
        }

        System.out.println("Requested " + hitProteins.size() + " proteins.");

        MessageStatus status = null;
        status = new MessageStatus("Sucess", 0, 0, "", "");
        return new MutablePair<>(result, status);
    }

    /*
     * This method expects the line to be validated already
     */
    private static Snp getSnpFromChrBp(String line) {
        String[] fields = line.split("\\s");
        Integer chr = Integer.valueOf(fields[0]);
        Long bp = Long.valueOf(fields[1]);

        return new Snp(chr, bp);
    }

    public static Pair<List<String[]>, MessageStatus> searchWithVCF(
            List<String> input,
            ImmutableMap<String, String> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, String> imapChrBpToProteins,
            ImmutableSetMultimap<String, String> imapProteinsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            HashSet<String> hitProteins,
            HashSet<String> hitPathways) {

        List<String[]> result = new ArrayList<String[]>();

        int row = 0;
        for (String line : input) {
            line = line.trim();
            row++;
            if (line.startsWith("#")) {
                continue;
            }
            if (line.isEmpty()) {
                sendWarning(EMPTY_ROW, row);
                continue;
            }
            if (!matches_Vcf_Record(line)) {
                sendWarning(INVALID_ROW, row);
                continue;
            }

            Snp snp = getSnpFromChrBp(line);
            for (String protein : imapChrBpToProteins.get(snp.getChr() + "_" + snp.getBp())) {
                for (String reaction : imapProteinsToReactions.get(protein)) {
                    for (String pathwayStId : imapReactionsToPathways.get(reaction)) {
                        if (topLevelPathways && imapPathwaysToTopLevelPathways.get(pathwayStId).size() > 0) {
                            for (String topLevelPathway : imapPathwaysToTopLevelPathways.get(pathwayStId)) {
                                String[] values = {
                                        snp.getChr() + "+" + snp.getBp(),
                                        protein,
                                        reaction,
                                        iReactions.get(reaction),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName(),
                                        topLevelPathway,
                                        iPathways.get(topLevelPathway).getDisplayName()
                                };
                                result.add(values);
                            }
                        } else {
                            String[] values = {
                                    snp.getChr() + "+" + snp.getBp(),
                                    protein,
                                    reaction,
                                    iReactions.get(reaction),
                                    pathwayStId,
                                    iPathways.get(pathwayStId).getDisplayName()
                            };
                            result.add(values);
                        }
                    }
                }
            }
        }

        System.out.println("Requested " + hitProteins.size() + " proteins.");

        MessageStatus status = null;
        status = new MessageStatus("Sucess", 0, 0, "", "");
        return new MutablePair<>(result, status);
    }

    public static Pair<List<String[]>, MessageStatus> searchWithProteoform(
            List<String> input,
            MatchType matchType,
            Long margin,
            ImmutableMap<String, String> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, Proteoform> imapProteinsToProteoforms,
            ImmutableSetMultimap<Proteoform, String> imapProteoformsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            HashSet<String> hitProteins,
            HashSet<String> hitPathways) {

        List<String[]> result = new ArrayList<String[]>();
        HashSet<Proteoform> inputProteoforms = new HashSet<>();
        HashSet<Proteoform> hitProteoforms = new HashSet<>();
        ProteoformMatcher matcher = null;
        switch (matchType) {
            case FLEXIBLE:
                matcher = new ProteoformMatcherFlexible();
                break;
            case ONE:
                matcher = new ProteoformMatcherOne();
                break;
            case STRICT:
                matcher = new ProteoformMatcherStrict();
                break;
        }

        int row = 1;
        for (String line : input) {
            row++;
            if (matches_Proteoform_Simple(line)) {
                try {
                    Proteoform proteoform = ProteoformFormat.SIMPLE.getProteoform(line, row);
                    inputProteoforms.add(proteoform);
                } catch (ParseException e) {
                    sendWarning(INVALID_ROW, row);
                }
            } else {
                if (line.isEmpty())
                    sendWarning(EMPTY_ROW, row);
                else
                    sendWarning(INVALID_ROW, row);
            }
        }

        for (Proteoform proteoform : inputProteoforms) {
            for (Proteoform refProteoform : imapProteinsToProteoforms.get(proteoform.getUniProtAcc())) {
                if (matcher.matches(proteoform, refProteoform, margin)) {
                    hitProteoforms.add(refProteoform);
                }
            }
        }

        for (Proteoform proteoform : hitProteoforms) {
            for (String reaction : imapProteoformsToReactions.get(proteoform)) {
                for (String pathway : imapReactionsToPathways.get(reaction)) {
                    if (topLevelPathways && imapPathwaysToTopLevelPathways.get(pathway).size() > 0) {
                        for (String topLevelPathway : imapPathwaysToTopLevelPathways.get(pathway)) {
                            String[] values = {
                                    proteoform.getUniProtAcc(),
                                    proteoform.toString(ProteoformFormat.SIMPLE),
                                    reaction,
                                    iReactions.get(reaction),
                                    pathway,
                                    iPathways.get(pathway).getDisplayName(),
                                    topLevelPathway,
                                    iPathways.get(topLevelPathway).getDisplayName()
                            };
                            result.add(values);
                        }
                    } else {
                        String[] values = {
                                proteoform.getUniProtAcc(),
                                proteoform.toString(ProteoformFormat.SIMPLE),
                                reaction,
                                iReactions.get(reaction),
                                pathway,
                                iPathways.get(pathway).getDisplayName()
                        };
                        result.add(values);
                    }
                }
            }
        }

        MessageStatus status = null;
        status = new MessageStatus("Sucess", 0, 0, "", "");
        return new MutablePair<>(result, status);
    }

    public static Pair<List<String[]>, MessageStatus> searchWithPeptide(
            Collection<String> input,
            ImmutableMap<String, String> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, String> imapProteinsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            HashSet<String> hitProteins,
            HashSet<String> hitPathways) {

        List<String[]> result = new ArrayList<String[]>();

        // Note: In this function the duplicate protein identifiers are removed by
        // adding the whole input list to a set.
        if (!initializePeptideMapper(fasta)) {
            System.out.println(ERROR_INITIALIZING_PEPTIDE_MAPPER.getMessage());
            System.exit(ERROR_INITIALIZING_PEPTIDE_MAPPER.getCode());
        }

        int row = 1;
        for (String line : input) {
            line = line.trim();
            row++;
            if (matches_Peptite(line)) {
                // Process line
                for (String protein : getPeptideMapping(line)) {
                    hitProteins.add(protein);
                }
            } else {
                if (line.isEmpty())
                    sendWarning(EMPTY_ROW, row);
                else
                    sendWarning(INVALID_ROW, row);
            }
        }

        return searchWithUniProt(hitProteins, iReactions, iPathways, imapProteinsToReactions, imapReactionsToPathways, imapPathwaysToTopLevelPathways, topLevelPathways, hitProteins, hitPathways);
    }

    public static Pair<List<String[]>, MessageStatus> searchWithModifiedPeptide(
            List<String> input,
            MatchType matchType,
            Long margin,
            ImmutableMap<String, String> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, Proteoform> imapProteinsToProteoforms,
            ImmutableSetMultimap<Proteoform, String> imapProteoformsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            HashSet<String> hitProteins,
            HashSet<String> hitPathways) {

        List<String[]> result = new ArrayList<String[]>();
        HashSet<Proteoform> inputProteoforms = new HashSet<>();
        HashSet<Proteoform> hitProteoforms = new HashSet<>();
        ProteoformMatcher matcher = null;

        // Note: In this function the duplicate protein identifiers are removed by
        // adding the whole input list to a set.
        if (!initializePeptideMapper(fasta)) {
            System.out.println(ERROR_INITIALIZING_PEPTIDE_MAPPER.getMessage());
            System.exit(ERROR_INITIALIZING_PEPTIDE_MAPPER.getCode());
        }

        int row = 1;
        for (String line : input) {
            row++;
            if (matches_Peptite_And_Mod_Sites(line)) {
                // Process line
                for (String protein : getPeptideMapping(line)) {
                    try {
                        inputProteoforms.addAll(getProteoforms(line));
                    } catch (ParseException e) {
                        e.printStackTrace(); //TODO Replace
                    }
                }
            } else {
                if (line.isEmpty())
                    sendWarning(EMPTY_ROW, row);
                else
                    sendWarning(INVALID_ROW, row);
            }
        }

        switch (matchType) {
            case FLEXIBLE:
                matcher = new ProteoformMatcherFlexible();
                break;
            case ONE:
                matcher = new ProteoformMatcherOne();
                break;
            case STRICT:
                matcher = new ProteoformMatcherStrict();
                break;
        }

        for (Proteoform proteoform : inputProteoforms) {
            for (Proteoform refProteoform : imapProteinsToProteoforms.get(proteoform.getUniProtAcc())) {
                if (matcher.matches(proteoform, refProteoform, margin)) {
                    hitProteoforms.add(refProteoform);
                }
            }
        }

        for (Proteoform proteoform : hitProteoforms) {
            for (String reaction : imapProteoformsToReactions.get(proteoform)) {
                for (String pathway : imapReactionsToPathways.get(reaction)) {
                    if (topLevelPathways && imapPathwaysToTopLevelPathways.get(pathway).size() > 0) {
                        for (String topLevelPathway : imapPathwaysToTopLevelPathways.get(pathway)) {
                            String[] values = {
                                    proteoform.getUniProtAcc(),
                                    proteoform.toString(ProteoformFormat.SIMPLE),
                                    reaction,
                                    iReactions.get(reaction),
                                    pathway,
                                    iPathways.get(pathway).getDisplayName(),
                                    topLevelPathway,
                                    iPathways.get(topLevelPathway).getDisplayName()
                            };
                            result.add(values);
                        }
                    } else {
                        String[] values = {
                                proteoform.getUniProtAcc(),
                                proteoform.toString(ProteoformFormat.SIMPLE),
                                reaction,
                                iReactions.get(reaction),
                                pathway,
                                iPathways.get(pathway).getDisplayName()
                        };
                        result.add(values);
                    }
                }
            }
        }

        MessageStatus status = null;
        status = new MessageStatus("Sucess", 0, 0, "", "");
        return new MutablePair<>(result, status);
    }
}
