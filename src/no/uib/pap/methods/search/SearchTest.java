package no.uib.pap.methods.search;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.io.Files;
import no.uib.pap.model.*;
import org.apache.commons.lang3.tuple.Pair;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.nio.charset.Charset;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.zip.GZIPInputStream;

import static org.junit.jupiter.api.Assertions.*;

class SearchTest {

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

    @Test
    void searchWithUniProtFillHitsTest1() {
        List<String> input = new ArrayList<>();
        input.add("P01308");
        input.add("P0130898989");

        Pair<List<String[]>, MessageStatus> result = Search.searchWithUniProt(
                input,
                iReactions,
                iPathways,
                imapProteinsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                true,
                hitProteins,
                hitPathways
        );

        assertEquals(1, hitProteins.size());
        assertTrue(hitProteins.contains("P01308"));
        assertEquals(29, hitPathways.size());
        assertTrue(hitPathways.contains("R-HSA-264876"));
        assertTrue(hitPathways.contains("R-HSA-74749"));
    }

    @Test
    void searchWithUniProtFillHitsTest2() {
        List<String> input = new ArrayList<>();
        input.add("P31749");
        input.add("blabla");
        input.add("Q5S007 ");
        input.add("P10636");

        Pair<List<String[]>, MessageStatus> result = Search.searchWithUniProt(
                input,
                iReactions,
                iPathways,
                imapProteinsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                true,
                hitProteins,
                hitPathways
        );

        assertEquals(3, hitProteins.size());
        assertTrue(hitProteins.contains("P10636"));
        assertEquals(92, hitPathways.size());
        assertTrue(hitPathways.contains("R-HSA-8857538"));
        assertTrue(hitPathways.contains("R-HSA-5663202"));
    }

    @Test
    void searchUniProtCountReactionsAndEntitiesFoundTest() throws ParseException {
        List<String> input = new ArrayList<>();
        input.add("P01308");

        Pair<List<String[]>, MessageStatus> result = Search.searchWithUniProt(
                input,
                iReactions,
                iPathways,
                imapProteinsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                true,
                hitProteins,
                hitPathways
        );

        assertEquals(1, iPathways.get("R-HSA-74749").getEntitiesFound().size());
        assertEquals(1, iPathways.get("R-HSA-74749").getReactionsFound().size());
        assertTrue(iPathways.get("R-HSA-74749").getReactionsFound().contains("R-HSA-110011"));

        assertEquals(6, iPathways.get("R-HSA-6807878").getReactionsFound().size());
        assertEquals(1, iPathways.get("R-HSA-6807878").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-6807878").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308")));
        assertTrue(iPathways.get("R-HSA-6807878").getReactionsFound().contains("R-HSA-6807875"));
        assertTrue(iPathways.get("R-HSA-6807878").getReactionsFound().contains("R-HSA-6809010"));
        assertTrue(iPathways.get("R-HSA-6807878").getReactionsFound().contains("R-HSA-6809003"));
    }

    @Test
    void searchWithGeneFillHitsTest1() {
        List<String> input = new ArrayList<>();
        input.add("GCK");
        input.add("GCK ");
        input.add("HNF4A");
        input.add("blabla");

        Pair<List<String[]>, MessageStatus> result = Search.searchWithGene(
                input,
                iReactions,
                iPathways,
                imapGenesToProteins,
                imapProteinsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                true,
                hitProteins,
                hitPathways
        );

        assertEquals(2, hitProteins.size());
        assertTrue(hitProteins.contains("P35557"));
        assertTrue(hitProteins.contains("P41235"));
        assertEquals(16, hitPathways.size());
        assertTrue(hitPathways.contains("R-HSA-170822"));
        assertTrue(hitPathways.contains("R-HSA-74160"));
        assertTrue(hitPathways.contains("R-HSA-1266738"));

    }

    @Test
    void searchWithGeneFillHitsTest2() throws ParseException {
        List<String> input = new ArrayList<>();
        input.add("CFTR");
        input.add("TGFB1 ");
        input.add("FCGR2A");
        input.add("DCTN4");
        input.add("SCNN1B");
        input.add("SCNN1G");
        input.add("SCNN1A");
        input.add("TNFRSF1A");
        input.add("CLCA4");
        input.add("STX1A");
        input.add("CXCL8");


        Pair<List<String[]>, MessageStatus> result = Search.searchWithGene(
                input,
                iReactions,
                iPathways,
                imapGenesToProteins,
                imapProteinsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                true,
                hitProteins,
                hitPathways
        );

        assertEquals(11, hitProteins.size());
        assertTrue(hitProteins.contains("P37088"));
        assertTrue(hitProteins.contains("P01137"));
        assertEquals(125, hitPathways.size());
        assertTrue(hitPathways.contains("R-HSA-2672351"));
        assertTrue(hitPathways.contains("R-HSA-76002"));
        assertTrue(hitPathways.contains("R-HSA-449147"));

        // Check counts
        assertEquals(4, iPathways.get("R-HSA-2672351").getEntitiesFound().size());
        assertEquals(6, iPathways.get("R-HSA-2672351").getReactionsFound().size());
        assertTrue(iPathways.get("R-HSA-2672351").getReactionsFound().contains("R-HSA-5333671"));
        assertTrue(iPathways.get("R-HSA-2672351").getReactionsFound().contains("R-HSA-2672334"));
        assertTrue(iPathways.get("R-HSA-2672351").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("Q14CN2")));
        assertTrue(iPathways.get("R-HSA-2672351").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P51170")));

        assertEquals(3, iPathways.get("R-HSA-1266738").getReactionsFound().size());
        assertEquals(2, iPathways.get("R-HSA-1266738").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-1266738").getReactionsFound().contains("R-HSA-381283"));
        assertTrue(iPathways.get("R-HSA-1266738").getReactionsFound().contains("R-HSA-560491"));
        assertTrue(iPathways.get("R-HSA-1266738").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("Q16623")));
        assertTrue(iPathways.get("R-HSA-1266738").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01137")));
    }

    @Test
    void searchWithEnsembl() {
    }

    @Test
    void searchWithRsId() {
    }

    @Test
    void searchWithChrBp() {
    }

    @Test
    void searchWithVCF() {
    }

    @Test
    void searchWithProteoform() throws IOException {
        Pair<List<String[]>, MessageStatus> result = Search.searchWithProteoform(
                Files.readLines(new File(resourcesPath + "input/ReactomeAllProteoformsSimple.csv"), Charset.defaultCharset()),
                MatchType.FLEXIBLE,
                0L,
                (ImmutableMap<String, String>) getSerializedObject("iReactions.gz"),
                (ImmutableMap<String, Pathway>) getSerializedObject("iPathways.gz"),
                (ImmutableSetMultimap<String, Proteoform>) getSerializedObject("imapProteinsToProteoforms.gz"),
                (ImmutableSetMultimap<Proteoform, String>) getSerializedObject("imapProteoformsToReactions.gz"),
                (ImmutableSetMultimap<String, String>) getSerializedObject("imapReactionsToPathways.gz"),
                (ImmutableSetMultimap<String, String>) getSerializedObject("imapPathwaysToTopLevelPathways.gz"),
                true,
                hitProteins,
                hitPathways
        );

        assertEquals(441275, result.getLeft().size());
    }

    @Test
    void searchWithPeptide() {
    }

    @Test
    void searchWithModifiedPeptide() {
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