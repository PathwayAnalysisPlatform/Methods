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
import java.util.HashSet;
import java.util.List;
import java.util.zip.GZIPInputStream;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

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
    void searchWithEnsemblFillHitsDiabetesInYouthTest() throws ParseException {
        List<String> input = new ArrayList<>();
        input.add("blabla");
        input.add("ENSG00000101076");
        input.add("ENSG00000106633");
        input.add("ENSP00000223366");
        input.add("ENSP00000312987");
        input.add("ENSP00000315180");
        input.add("ENSP00000379142");
        input.add("ENSP00000384247");
        input.add("ENSP00000396216");
        input.add("ENSP00000410911");
        input.add("\t\tENSP00000412111");
        input.add("ENSP00000476609");
        input.add("ENSP00000482149   ");


        Pair<List<String[]>, MessageStatus> result = Search.searchWithEnsembl(
                input,
                iReactions,
                iPathways,
                imapEnsemblToProteins,
                imapProteinsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                true,
                hitProteins,
                hitPathways
        );

        assertEquals(2, hitProteins.size());
        assertTrue(hitProteins.contains("P41235"));
        assertTrue(hitProteins.contains("P35557"));
        assertEquals(16, hitPathways.size());
        assertTrue(hitPathways.contains("R-HSA-170822"));
        assertTrue(hitPathways.contains("R-HSA-74160"));
        assertTrue(hitPathways.contains("R-HSA-1266738"));

        // Check counts
        assertEquals(5, iPathways.get("R-HSA-170822").getReactionsFound().size());
        assertTrue(iPathways.get("R-HSA-170822").getReactionsFound().contains("R-HSA-170796"));
        assertTrue(iPathways.get("R-HSA-170822").getReactionsFound().contains("R-HSA-170810"));
        assertTrue(iPathways.get("R-HSA-170822").getReactionsFound().contains("R-HSA-170825"));
        assertEquals(1, iPathways.get("R-HSA-170822").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-170822").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P35557")));

        assertEquals(1, iPathways.get("R-HSA-383280").getReactionsFound().size());
        assertTrue(iPathways.get("R-HSA-383280").getReactionsFound().contains("R-HSA-376419"));
        assertEquals(1, iPathways.get("R-HSA-383280").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-383280").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P41235")));
    }

    @Test
    void searchWithEnsemblFillHitsCysticFibrosisTest() throws ParseException {
        List<String> input = new ArrayList<>();
        input.add("ENSG00000001626");
        input.add("ENSG00000016602");
        input.add("ENSG00000067182");
        input.add("ENSG00000105329");
        input.add("ENSG00000106089");
        input.add("\tENSG00000111319");
        input.add("ENSG00000132912");
        input.add("ENSG00000143226");
        input.add("   ENSG00000166828");
        input.add("ENSG00000168447");
        input.add("ENSG00000169429");
        input.add("blabla");

        Pair<List<String[]>, MessageStatus> result = Search.searchWithEnsembl(
                input,
                iReactions,
                iPathways,
                imapEnsemblToProteins,
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
    void searchWithRsId() {
    }

    @Test
    void searchWithChrBp() {
    }

    @Test
    void searchWithVCF() {
    }

    @Test
    void searchWithProteoformTest() throws IOException {
        Pair<List<String[]>, MessageStatus> result = Search.searchWithProteoform(
                Files.readLines(new File(resourcesPath + "input/ReactomeAllProteoformsSimple.csv"), Charset.defaultCharset()),
                MatchType.FLEXIBLE,
                0L,
                iReactions,
                iPathways,
                imapProteinsToProteoforms,
                imapProteoformsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                true,
                hitProteins,
                hitPathways
        );

        assertEquals(441275, result.getLeft().size());
    }

    @Test
    // All proteoforms of a protein
    void searchWithProteoformSet1FillHitsTest() throws IOException {
        Pair<List<String[]>, MessageStatus> result = Search.searchWithProteoform(
                Files.readLines(new File(resourcesPath + "input/Proteoforms/Set1.csv"), Charset.defaultCharset()),
                MatchType.FLEXIBLE,
                0L,
                iReactions,
                iPathways,
                imapProteinsToProteoforms,
                imapProteoformsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                true,
                hitProteins,
                hitPathways
        );

//        assert()
    }

    @Test
    // Some proteoforms of a protein
    void searchWithProteoformSet2FillHitsTest() throws IOException {
        Pair<List<String[]>, MessageStatus> result = Search.searchWithProteoform(
                Files.readLines(new File(resourcesPath + "input/Proteoforms/Set2.csv"), Charset.defaultCharset()),
                MatchType.FLEXIBLE,
                0L,
                iReactions,
                iPathways,
                imapProteinsToProteoforms,
                imapProteoformsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                true,
                hitProteins,
                hitPathways
        );

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