package no.uib.pap.methods.analysis.ora;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSetMultimap;
import no.uib.pap.methods.search.Search;
import no.uib.pap.model.MessageStatus;
import no.uib.pap.model.Pathway;
import no.uib.pap.model.Proteoform;
import no.uib.pap.model.ProteoformFormat;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;

import java.io.InputStream;
import java.io.ObjectInputStream;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

import static org.junit.jupiter.api.Assertions.*;

class AnalysisTest {

    static String resourcesPath = "../../PathwayMatcher/resources/";
    static HashSet<String> hitProteins = null;
    static HashSet<String> hitPathways = null;

    private static ImmutableMap<String, String> iReactions; // Reaction stId to Reaction displayName
    private static ImmutableMap<String, Pathway> iPathways; // Pathway stId to Pathway instance
    private static ImmutableMap<String, String> iProteins; // Protein accession (UniProt) to name
    private static ImmutableSetMultimap<String, String> imapGenesToProteins = null;
    private static ImmutableSetMultimap<String, String> imapEnsemblToProteins = null;
    private static ImmutableSetMultimap<String, String> imapProteinsToReactions = null;
    private static ImmutableSetMultimap<String, String> imapReactionsToPathways = null;
    private static ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways = null;
    private static ImmutableSetMultimap<String, Proteoform> imapProteinsToProteoforms = null;
    private static ImmutableSetMultimap<Proteoform, String> imapProteoformsToReactions = null;
    private static ImmutableSetMultimap<String, String> imapRsIdsToProteins = null;
    private static ImmutableSetMultimap<String, String> imapChrBpToProteins = null;
    private static ImmutableSetMultimap<String, String> imapReactionsToParticipants = null;
    private static ImmutableSetMultimap<String, String> imapProteinsToComplexes = null;
    private static ImmutableSetMultimap<String, String> imapComplexesToParticipants = null;

    @BeforeAll
    static void loadStaticMapping() {
        iReactions = (ImmutableMap<String, String>) getSerializedObject("iReactions.gz");
        iPathways = (ImmutableMap<String, Pathway>) getSerializedObject("iPathways.gz");
        iProteins = (ImmutableMap<String, String>) getSerializedObject("iProteins.gz");
        imapGenesToProteins = (ImmutableSetMultimap<String, String>) getSerializedObject("imapGenesToProteins.gz");
        imapEnsemblToProteins = (ImmutableSetMultimap<String, String>) getSerializedObject("imapEnsemblToProteins.gz");
        imapProteinsToReactions = (ImmutableSetMultimap<String, String>) getSerializedObject("imapProteinsToReactions.gz");
        imapReactionsToPathways = (ImmutableSetMultimap<String, String>) getSerializedObject("imapReactionsToPathways.gz");
        imapPathwaysToTopLevelPathways = (ImmutableSetMultimap<String, String>) getSerializedObject("imapPathwaysToTopLevelPathways.gz");
        imapProteinsToProteoforms = (ImmutableSetMultimap<String, Proteoform>) getSerializedObject("imapProteinsToProteoforms.gz");
        imapProteoformsToReactions = (ImmutableSetMultimap<Proteoform, String>) getSerializedObject("imapProteoformsToReactions.gz");
//        imapRsIdsToProteins = (ImmutableSetMultimap<String,String>) getSerializedObject("imapRsIdsToProteins.gz");
//        imapChrBpToProteins = (ImmutableSetMultimap<String,String>) getSerializedObject("imapChrBpToProteins.gz");
        imapReactionsToParticipants = (ImmutableSetMultimap<String, String>) getSerializedObject("imapReactionsToParticipants.gz");
        imapProteinsToComplexes = (ImmutableSetMultimap<String, String>) getSerializedObject("imapProteinsToComplexes.gz");
        imapComplexesToParticipants = (ImmutableSetMultimap<String, String>) getSerializedObject("imapComplexesToParticipants.gz");
    }

    @BeforeEach
    void setUp() {
        hitProteins = new HashSet<>();
        hitPathways = new HashSet<>();
    }

    @AfterEach
    void tearDown() {
        hitProteins.clear();
        hitPathways.clear();
    }

    @org.junit.jupiter.api.Test
    void analysisTest() throws ParseException {
        MessageStatus messageStatus = new MessageStatus("Sucess", 0, 0, "", "");
        HashSet<String> hitPathways = new HashSet<>();
        HashSet<String> hitProteins = new HashSet<>();
        List<String> input = new ArrayList<>();

        input.add("P01308");
        // Execute the search to fill the iPathways, hitPahtways and hitProteins data structures
        Search.searchWithUniProt(input,
                iReactions,
                iPathways,
                imapProteinsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                true,
                hitProteins,
                hitPathways);

        messageStatus = Analysis.analysis(iPathways, imapProteinsToReactions.keySet().size(), hitProteins, hitPathways);

        assertEquals(29, hitPathways.size());
        assertTrue(hitPathways.contains("R-HSA-1430728"));
        assertEquals(1, iPathways.get("R-HSA-74749").getReactionsFound().size());
        assertEquals(1, iPathways.get("R-HSA-74749").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-74749").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308")));
        assertTrue(iPathways.get("R-HSA-74749").getReactionsFound().contains("R-HSA-110011"));

        assertEquals(0.001349, iPathways.get("R-HSA-5653656").getReactionsRatio(), 0.01);
        assertEquals(0.001199, iPathways.get("R-HSA-5653656").getEntitiesRatio(), 0.01);
    }

    public static Object getSerializedObject(String fileName) {
        Object obj = null;
        try {
            InputStream fileStream = ClassLoader.getSystemResourceAsStream(fileName);
            GZIPInputStream gis = new GZIPInputStream(fileStream);
            ObjectInputStream ois = new ObjectInputStream(gis);
            obj = ois.readObject();
            ois.close();

        } catch (Exception ex) {
            System.out.println("Error loading file: " + fileName);
            ex.printStackTrace();
        }
        return obj;
    }
}