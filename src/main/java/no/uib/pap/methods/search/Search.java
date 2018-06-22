package no.uib.pap.methods.search;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSetMultimap;
import no.uib.pap.model.*;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;

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

    // Fills the hitProteins set to call the next method
    public static Pair<List<String[]>, MessageStatus> searchWithUniProt(
            Collection<String> input,
            ImmutableMap<String, Reaction> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, String> imapProteinsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            TreeSet<String> hitProteins,
            HashSet<String> hitPathways) {

        List<String[]> result = new ArrayList<String[]>();

        for (String protein : input) {
            protein = protein.trim();

            if (protein.contains("-")) {
                protein = protein.substring(0, protein.indexOf("-"));
            }

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
                    if (topLevelPathways) {
                        if (imapPathwaysToTopLevelPathways.get(pathwayStId).size() > 0) {
                            for (String topLevelPathway : imapPathwaysToTopLevelPathways.get(pathwayStId)) {
                                String[] values = {
                                        protein,
                                        reaction,
                                        iReactions.get(reaction).getDisplayName(),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName(),
                                        topLevelPathway,
                                        iPathways.get(topLevelPathway).getDisplayName()
                                };
                                result.add(values);
                            }
                        } else {
                            String[] values = {
                                    protein,
                                    reaction,
                                    iReactions.get(reaction).getDisplayName(),
                                    pathwayStId,
                                    iPathways.get(pathwayStId).getDisplayName(),
                                    pathwayStId,
                                    iPathways.get(pathwayStId).getDisplayName()
                            };
                            result.add(values);
                        }
                    } else {
                        String[] values = {
                                protein,
                                reaction,
                                iReactions.get(reaction).getDisplayName(),
                                pathwayStId,
                                iPathways.get(pathwayStId).getDisplayName()
                        };
                        result.add(values);
                    }
                }
            }
        }

        System.out.println("\nRequested " + hitProteins.size() + " proteins.");

        MessageStatus status = null;
        status = new MessageStatus("Success", 0, 0, "", "");
        return new MutablePair<>(result, status);
    }


    public static Pair<List<String[]>, MessageStatus> searchWithGene(
            List<String> input,
            ImmutableMap<String, Reaction> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, String> imapGenesToProteins,
            ImmutableSetMultimap<String, String> imapProteinsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            TreeSet<String> hitProteins,
            HashSet<String> hitPathways,
            TreeSet<String> hitGenes) {

        List<String[]> result = new ArrayList<String[]>();

        for (String gene : input) {
            for (String protein : imapGenesToProteins.get(gene.trim())) {
                hitGenes.add(gene);
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
                        if (topLevelPathways) {
                            if (imapPathwaysToTopLevelPathways.get(pathwayStId).size() > 0) {
                                for (String topLevelPathway : imapPathwaysToTopLevelPathways.get(pathwayStId)) {
                                    String[] values = {
                                            gene,
                                            protein,
                                            reaction,
                                            iReactions.get(reaction).getDisplayName(),
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
                                        iReactions.get(reaction).getDisplayName(),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName(),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName()
                                };
                                result.add(values);
                            }
                        } else {
                            String[] values = {
                                    gene,
                                    protein,
                                    reaction,
                                    iReactions.get(reaction).getDisplayName(),
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
        status = new MessageStatus("Success", 0, 0, "", "");
        return new MutablePair<>(result, status);
    }

    public static Pair<List<String[]>, MessageStatus> searchWithEnsembl(
            List<String> input,
            ImmutableMap<String, Reaction> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, String> imapEnsemblToProteins,
            ImmutableSetMultimap<String, String> imapProteinsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            TreeSet<String> hitProteins,
            HashSet<String> hitPathways) {

        List<String[]> result = new ArrayList<String[]>();

        for (String ensembl : input) {
            for (String protein : imapEnsemblToProteins.get(ensembl.trim())) {
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
                        if (topLevelPathways) {
                            if (imapPathwaysToTopLevelPathways.get(pathwayStId).size() > 0) {
                                for (String topLevelPathway : imapPathwaysToTopLevelPathways.get(pathwayStId)) {
                                    String[] values = {
                                            ensembl,
                                            protein,
                                            reaction,
                                            iReactions.get(reaction).getDisplayName(),
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
                                        iReactions.get(reaction).getDisplayName(),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName(),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName(),
                                };
                                result.add(values);
                            }
                        } else {
                            String[] values = {
                                    ensembl,
                                    protein,
                                    reaction,
                                    iReactions.get(reaction).getDisplayName(),
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
        status = new MessageStatus("Success", 0, 0, "", "");
        return new MutablePair<>(result, status);
    }


    /**
     * Maps rsids to protein to reaction to pathways.
     * Usually only for rsids in a specific chromosome, but uses all the mapping contained at the imapRsIdsToProteins parameter.
     * Fills the hitPathways and the hitProteins from the parameter structures.
     *
     * @param input                          Set of unique identifiers
     * @param iReactions                     Structure to read reactions stId and displayName
     * @param iPathways                      Here add the reactions and entities found
     * @param imapRsIdsToProteins            The mapping for one chromosome
     * @param imapProteinsToReactions        Generic mapping all proteins to all reactions
     * @param imapReactionsToPathways        Generic mapping all reactions to all pathways
     * @param imapPathwaysToTopLevelPathways Generic mapping all pathways to top level pathways
     * @param topLevelPathways               Flag if top level pathways should be used
     * @param hitProteins                    Structure to keep which proteins were hit by the search
     * @param hitPathways                    Structure to keep which pathways were hit by the search
     * @return Mapping from rsids to pathways, message errors
     */
    public static Pair<List<String[]>, MessageStatus> searchWithRsId(
            HashSet<String> input,
            ImmutableMap<String, Reaction> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, String> imapRsIdsToProteins,
            ImmutableSetMultimap<String, String> imapProteinsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            TreeSet<String> hitProteins,
            HashSet<String> hitPathways) {

        List<String[]> result = new ArrayList<String[]>();

        // If a variant is found then discard it from the variant set
        for (String rsid : input) {
//            if (imapRsIdsToProteins.containsKey(rsid)) {
//                System.out.println("Analysing: " + rsid);
//            }
            for (String protein : imapRsIdsToProteins.get(rsid)) {
                //System.out.println("Mapped to: " + protein);
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

                        if (topLevelPathways) {
                            if (imapPathwaysToTopLevelPathways.get(pathwayStId).size() > 0) {
                                for (String topLevelPathway : imapPathwaysToTopLevelPathways.get(pathwayStId)) {
                                    String[] values = {
                                            rsid,
                                            protein,
                                            reaction,
                                            iReactions.get(reaction).getDisplayName(),
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
                                        iReactions.get(reaction).getDisplayName(),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName(),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName(),
                                };
                                result.add(values);
                            }
                        } else {
                            String[] values = {
                                    rsid,
                                    protein,
                                    reaction,
                                    iReactions.get(reaction).getDisplayName(),
                                    pathwayStId,
                                    iPathways.get(pathwayStId).getDisplayName()
                            };
                            result.add(values);
                        }
                    }
                }
            }

        }

        System.out.println("Found " + hitProteins.size() + " proteins.");

        MessageStatus status = null;
        status = new MessageStatus("Success", 0, 0, "", "");
        return new MutablePair<>(result, status);
    }

    /**
     * Maps variants composed by [chr, bp] to protein to reaction to pathways.
     * Only maps a specific chromosome at a time, but uses all the mapping contained at the imapRsIdsToProteins parameter.
     * Fills the hitPathways and the hitProteins from the parameter structures.
     *
     * @param bpSet                          Set of base pairs in the desired chromosome
     * @param iReactions                     Structure to read reactions stId and displayName
     * @param iPathways                      Here add the reactions and entities found
     * @param imapChrBpToProteins            The mapping for one chromosome
     * @param imapProteinsToReactions        Generic mapping all proteins to all reactions
     * @param imapReactionsToPathways        Generic mapping all reactions to all pathways
     * @param imapPathwaysToTopLevelPathways Generic mapping all pathways to top level pathways
     * @param topLevelPathways               Flag if top level pathways should be used
     * @param hitProteins                    Structure to keep which proteins were hit by the search
     * @param hitPathways                    Structure to keep which pathways were hit by the search
     * @return Mapping from rsids to pathways, message errors
     */
    public static Pair<List<String[]>, MessageStatus> searchWithChrBp(
            int chr,
            Set<Long> bpSet,
            ImmutableMap<String, Reaction> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<Long, String> imapChrBpToProteins,
            ImmutableSetMultimap<String, String> imapProteinsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            TreeSet<String> hitProteins,
            HashSet<String> hitPathways) {

        List<String[]> result = new ArrayList<String[]>();

        int row = 0;
        for (Long bp : bpSet) {
            for (String protein : imapChrBpToProteins.get(bp)) {
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

                        if (topLevelPathways) {
                            if (imapPathwaysToTopLevelPathways.get(pathwayStId).size() > 0) {
                                for (String topLevelPathway : imapPathwaysToTopLevelPathways.get(pathwayStId)) {
                                    String[] values = {
                                            String.valueOf(chr),
                                            String.valueOf(bp),
                                            protein,
                                            reaction,
                                            iReactions.get(reaction).getDisplayName(),
                                            pathwayStId,
                                            iPathways.get(pathwayStId).getDisplayName(),
                                            topLevelPathway,
                                            iPathways.get(topLevelPathway).getDisplayName()
                                    };
                                    result.add(values);
                                }
                            } else {
                                String[] values = {
                                        String.valueOf(chr),
                                        String.valueOf(bp),
                                        protein,
                                        reaction,
                                        iReactions.get(reaction).getDisplayName(),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName(),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName(),
                                };
                                result.add(values);
                            }
                        } else {
                            String[] values = {
                                    String.valueOf(chr),
                                    String.valueOf(bp),
                                    protein,
                                    reaction,
                                    iReactions.get(reaction).getDisplayName(),
                                    pathwayStId,
                                    iPathways.get(pathwayStId).getDisplayName()
                            };
                            result.add(values);
                        }
                    }
                }
            }
        }

        System.out.println("Found " + hitProteins.size() + " proteins.");

        MessageStatus status = null;
        status = new MessageStatus("Success", 0, 0, "", "");
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

    public static Pair<List<String[]>, MessageStatus> searchWithProteoform(
            List<String> input,
            MatchType matchType,
            Long margin,
            ImmutableMap<String, Reaction> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, Proteoform> imapProteinsToProteoforms,
            ImmutableSetMultimap<Proteoform, String> imapProteoformsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            TreeSet<String> hitProteins,
            HashSet<Proteoform> hitProteoforms,
            HashSet<String> hitPathways) {

        List<String[]> result = new ArrayList<String[]>();
        HashSet<Proteoform> inputProteoforms = new HashSet<>();
        ProteoformMatching matcher = ProteoformMatching.getInstance(matchType);

        assert matcher != null;

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
                hitProteins.add(proteoform.getUniProtAcc());
                for (String pathwayStId : imapReactionsToPathways.get(reaction)) {
                    hitPathways.add(pathwayStId);
                    Pathway pathway = iPathways.get(pathwayStId);
                    pathway.getReactionsFound().add(reaction);
                    pathway.getEntitiesFound().add(proteoform);

                    if (topLevelPathways) {
                        if (imapPathwaysToTopLevelPathways.get(pathwayStId).size() > 0) {
                            for (String topLevelPathway : imapPathwaysToTopLevelPathways.get(pathwayStId)) {
                                String[] values = {
                                        proteoform.getUniProtAcc(),
                                        proteoform.toString(ProteoformFormat.SIMPLE),
                                        reaction,
                                        iReactions.get(reaction).getDisplayName(),
                                        pathwayStId,
                                        iPathways.get(pathwayStId).getDisplayName(),
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
                                    iReactions.get(reaction).getDisplayName(),
                                    pathwayStId,
                                    iPathways.get(pathwayStId).getDisplayName(),
                                    pathwayStId,
                                    iPathways.get(pathwayStId).getDisplayName(),
                            };
                            result.add(values);
                        }
                    } else {
                        String[] values = {
                                proteoform.getUniProtAcc(),
                                proteoform.toString(ProteoformFormat.SIMPLE),
                                reaction,
                                iReactions.get(reaction).getDisplayName(),
                                pathwayStId,
                                iPathways.get(pathwayStId).getDisplayName()
                        };
                        result.add(values);
                    }
                }
            }
        }

        MessageStatus status = null;
        status = new MessageStatus("Success", 0, 0, "", "");
        return new MutablePair<>(result, status);
    }

    public static Pair<List<String[]>, MessageStatus> searchWithPeptide(
            Collection<String> input,
            ImmutableMap<String, Reaction> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, String> imapProteinsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            TreeSet<String> hitProteins,
            HashSet<String> hitPathways,
            String fastaFile) {

        List<String[]> result = new ArrayList<String[]>();

        // Note: In this function the duplicate protein identifiers are removed by
        // adding the whole input list to a set.
        if (!initializePeptideMapper(fastaFile)) {
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
                    if (protein.contains("-")) {
                        hitProteins.add(protein.substring(0, protein.indexOf("-")));
                    } else {
                        hitProteins.add(protein);
                    }
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
            ImmutableMap<String, Reaction> iReactions,
            ImmutableMap<String, Pathway> iPathways,
            ImmutableSetMultimap<String, Proteoform> imapProteinsToProteoforms,
            ImmutableSetMultimap<Proteoform, String> imapProteoformsToReactions,
            ImmutableSetMultimap<String, String> imapReactionsToPathways,
            ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways,
            Boolean topLevelPathways,
            TreeSet<String> hitProteins,
            HashSet<Proteoform> hitProteoforms,
            HashSet<String> hitPathways,
            String fastaFile) {

        List<String[]> result = new ArrayList<String[]>();
        HashSet<Proteoform> inputProteoforms = new HashSet<>();
        ProteoformMatching matcher = ProteoformMatching.getInstance(matchType);
        List<String> correctedInput = new ArrayList<>();

        // Note: In this function the duplicate protein identifiers are removed by
        // adding the whole input list to a set.
        if (!initializePeptideMapper(fastaFile)) {
            System.out.println(ERROR_INITIALIZING_PEPTIDE_MAPPER.getMessage());
            System.exit(ERROR_INITIALIZING_PEPTIDE_MAPPER.getCode());
        }

        int row = 1;
        for (String line : input) {
            row++;
            if (matches_Peptite_And_Mod_Sites(line)) {
                try {
                    Proteoform tempProteoform = ProteoformFormat.SIMPLE.getProteoform(line);
                    for (Pair<String, Integer> pair : getPeptideMappingWithIndex(tempProteoform.getUniProtAcc())) {

                        String uniprot = pair.getLeft();
                        int index = pair.getRight();
                        Proteoform correctProteoform = new Proteoform(uniprot);

                        //Correct the positions of the PTMs
                        for(Pair<String, Long> ptm : tempProteoform.getPtms()){
                            correctProteoform.addPtm(ptm.getLeft(), ptm.getValue()+index);
                        }

                        correctedInput.add(correctProteoform.toString(ProteoformFormat.SIMPLE));
                    }
                } catch (ParseException e) {
                    e.printStackTrace();
                }
            } else {
                if (line.isEmpty())
                    sendWarning(EMPTY_ROW, row);
                else
                    sendWarning(INVALID_ROW, row);
            }
        }

        return searchWithProteoform(
                correctedInput,
                matchType,
                margin,
                iReactions,
                iPathways,
                imapProteinsToProteoforms,
                imapProteoformsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                topLevelPathways,
                hitProteins,
                hitProteoforms,
                hitPathways
        );
    }
}
