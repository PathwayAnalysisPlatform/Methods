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
import java.util.*;
import java.util.zip.GZIPInputStream;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class SearchTest {

    static String resourcesPath = "../../PathwayMatcher/resources/";
    static TreeSet<String> hitProteins = null;
    static HashSet<String> hitPathways = null;

    private static ImmutableMap<String, Reaction> iReactions; // Reaction stId to Reaction displayName
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
    private static ImmutableSetMultimap<Long, String> imapChrBpToProteins = null;
    private static ImmutableSetMultimap<String, String> imapReactionsToParticipants = null;

    @BeforeAll
    static void loadStaticMapping() {
        iReactions = (ImmutableMap<String, Reaction>) getSerializedObject("iReactions.gz");
        iPathways = (ImmutableMap<String, Pathway>) getSerializedObject("iPathways.gz");
        iProteins = (ImmutableMap<String, String>) getSerializedObject("iProteins.gz");
        imapGenesToProteins = (ImmutableSetMultimap<String, String>) getSerializedObject("imapGenesToProteins.gz");
        imapEnsemblToProteins = (ImmutableSetMultimap<String, String>) getSerializedObject("imapEnsemblToProteins.gz");
        imapProteinsToReactions = (ImmutableSetMultimap<String, String>) getSerializedObject("imapProteinsToReactions.gz");
        imapReactionsToPathways = (ImmutableSetMultimap<String, String>) getSerializedObject("imapReactionsToPathways.gz");
        imapPathwaysToTopLevelPathways = (ImmutableSetMultimap<String, String>) getSerializedObject("imapPathwaysToTopLevelPathways.gz");
        imapProteinsToProteoforms = (ImmutableSetMultimap<String, Proteoform>) getSerializedObject("imapProteinsToProteoforms.gz");
        imapProteoformsToReactions = (ImmutableSetMultimap<Proteoform, String>) getSerializedObject("imapProteoformsToReactions.gz");
    }

    @BeforeEach
    void setUp() {
        hitProteins = new TreeSet<>();
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
    void searchWithProteoformTest() throws IOException {
        Pair<List<String[]>, MessageStatus> result = Search.searchWithProteoform(
                Files.readLines(new File(resourcesPath + "input/ReactomeAllProteoformsSimple.csv"), Charset.defaultCharset()),
                MatchType.SUPERSET,
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

        assertEquals(456120, result.getLeft().size());
    }

    @Test
        // All proteoforms of a protein
    void searchWithProteoformSet1FillHitsTest() throws IOException, ParseException {
        Pair<List<String[]>, MessageStatus> result = Search.searchWithProteoform(
                Files.readLines(new File(resourcesPath + "input/Proteoforms/SIMPLE/Set1.csv"), Charset.defaultCharset()),
                MatchType.SUPERSET,
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

        assertEquals(1, hitProteins.size());
        assertEquals(29, hitPathways.size());

        assertEquals(3, iPathways.get("R-HSA-977225").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-977225").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00798:31,00798:43")));
        assertTrue(iPathways.get("R-HSA-977225").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00087:53,00798:31,00798:43")));
        assertFalse(iPathways.get("R-HSA-977225").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308")));

        assertEquals(1, iPathways.get("R-HSA-977225").getReactionsFound().size());
        assertTrue(iPathways.get("R-HSA-977225").getReactionsFound().contains("R-HSA-977136"));

        assertEquals(1, iPathways.get("R-HSA-199991").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-199991").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00798:31,00798:43,00798:95,00798:96,00798:100,00798:109")));
        assertFalse(iPathways.get("R-HSA-199991").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00087:53,00798:31,00798:43")));
        assertEquals(6, iPathways.get("R-HSA-199991").getReactionsFound().size());
        assertTrue(iPathways.get("R-HSA-199991").getReactionsFound().contains("R-HSA-6807877"));
        assertTrue(iPathways.get("R-HSA-199991").getReactionsFound().contains("R-HSA-6809003"));
    }

    @Test
    void searchWithProteoformsInsulinTest() throws IOException, ParseException {
        Pair<List<String[]>, MessageStatus> result = Search.searchWithProteoform(
                Files.readLines(new File(resourcesPath + "input/Proteoforms/SIMPLE/Insulin.txt"), Charset.defaultCharset()),
                MatchType.SUPERSET,
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

        assertEquals(1, hitProteins.size());
        assertEquals(22, hitPathways.size());

        assertEquals(3, iPathways.get("R-HSA-977225").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-977225").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00798:31,00798:43")));
        assertTrue(iPathways.get("R-HSA-977225").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00087:53,00798:31,00798:43")));
        assertFalse(iPathways.get("R-HSA-977225").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308")));


    }

    @Test
        // Some proteoforms of a protein
    void searchWithProteoformSet2FillHitsTest() throws IOException, ParseException {
        Pair<List<String[]>, MessageStatus> result = Search.searchWithProteoform(
                Files.readLines(new File(resourcesPath + "input/Proteoforms/SIMPLE/Set2.csv"), Charset.defaultCharset()),
                MatchType.SUPERSET,
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

        assertEquals(22, hitPathways.size());
        assertEquals(1, hitProteins.size());

        assertEquals(2, iPathways.get("R-HSA-977225").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-977225").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00798:31,00798:43")));
        assertTrue(iPathways.get("R-HSA-977225").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00087:53,00798:31,00798:43")));
        assertFalse(iPathways.get("R-HSA-977225").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308")));
        assertFalse(iPathways.get("R-HSA-977225").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00798:95,00798:96,00798:100,00798:109")));

        assertEquals(1, iPathways.get("R-HSA-977225").getReactionsFound().size());
        assertTrue(iPathways.get("R-HSA-977225").getReactionsFound().contains("R-HSA-977136"));

        assertEquals(3, iPathways.get("R-HSA-392499").getEntitiesFound().size());
        assertFalse(iPathways.get("R-HSA-977225").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00798:95,00798:96,00798:100,00798:109")));
        assertTrue(iPathways.get("R-HSA-392499").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00087:53,00798:31,00798:43")));
        assertTrue(iPathways.get("R-HSA-392499").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308")));

        assertEquals(10, iPathways.get("R-HSA-392499").getReactionsFound().size());
        assertTrue(iPathways.get("R-HSA-392499").getReactionsFound().contains("R-HSA-977136"));
        assertTrue(iPathways.get("R-HSA-392499").getReactionsFound().contains("R-HSA-9023178"));
    }

    @Test
    void singleProteoformSearchTest() throws IOException, ParseException {
        Pair<List<String[]>, MessageStatus> result = Search.searchWithProteoform(
                Files.readLines(new File(resourcesPath + "input/Proteoforms/SIMPLE/SingleProteoform.txt"), Charset.defaultCharset()),
                MatchType.SUPERSET,
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

        assertEquals(1, hitProteins.size());
        assertEquals(11, hitPathways.size());

        assertEquals(2, iPathways.get("R-HSA-168249").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;00048:127,00048:132,00048:171,00048:191,00048:226")));
        assertTrue(iPathways.get("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;")));
        assertFalse(iPathways.get("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561;00048:156,00048:161,00048:200,00048:220,00048:255")));
        assertFalse(iPathways.get("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561;00048:200,00048:220")));
    }

    @Test
    void singleProteoformSearchStrictTest() throws IOException, ParseException {
        Pair<List<String[]>, MessageStatus> result = Search.searchWithProteoform(
                Files.readLines(new File(resourcesPath + "input/Proteoforms/SIMPLE/SingleProteoform.txt"), Charset.defaultCharset()),
                MatchType.STRICT,
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

        assertEquals(1, hitProteins.size());
        assertEquals(11, hitPathways.size());

        assertEquals(1, iPathways.get("R-HSA-168249").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;00048:127,00048:132,00048:171,00048:191,00048:226")));
        assertFalse(iPathways.get("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;")));
    }

    @Test
    void singleModifiedPeptideSearchTest() throws IOException, ParseException {
        Pair<List<String[]>, MessageStatus> result = Search.searchWithModifiedPeptide(
                Files.readLines(new File(resourcesPath + "input/ModifiedPeptides/SingleModifiedPeptide.txt"), Charset.defaultCharset()),
                MatchType.SUPERSET,
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

        assertEquals(1, hitProteins.size());
        assertEquals(11, hitPathways.size());

        assertEquals(2, iPathways.get("R-HSA-168249").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;00048:127,00048:132,00048:171,00048:191,00048:226")));
        assertTrue(iPathways.get("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;")));
        assertFalse(iPathways.get("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561;00048:156,00048:161,00048:200,00048:220,00048:255")));
        assertFalse(iPathways.get("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561;00048:200,00048:220")));
    }

    @Test
    void singleModifiedPeptideSearchStrictTest() throws IOException, ParseException {
        Pair<List<String[]>, MessageStatus> result = Search.searchWithModifiedPeptide(
                Files.readLines(new File(resourcesPath + "input/ModifiedPeptides/SingleModifiedPeptide.txt"), Charset.defaultCharset()),
                MatchType.STRICT,
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

        assertEquals(1, hitProteins.size());
        assertEquals(11, hitPathways.size());

        assertEquals(1, iPathways.get("R-HSA-168249").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;00048:127,00048:132,00048:171,00048:191,00048:226")));
        assertFalse(iPathways.get("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;")));
    }

    @Test
    void singleModifiedPeptideDisplacedSearchTest() throws IOException, ParseException {
        Pair<List<String[]>, MessageStatus> result = Search.searchWithModifiedPeptide(
                Files.readLines(new File(resourcesPath + "input/ModifiedPeptides/SingleModifiedPeptideDisplaced.txt"), Charset.defaultCharset()),
                MatchType.STRICT,
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

        assertEquals(1, hitProteins.size());
        assertEquals(11, hitPathways.size());

        assertEquals(1, iPathways.get("R-HSA-168249").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;00048:127,00048:132,00048:171,00048:191,00048:226")));
        assertFalse(iPathways.get("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;")));
    }

    @Test
    void searchWithPeptideFillHitsTest1() {
        List<String> input = new ArrayList<>();
        input.add("LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN");

        Pair<List<String[]>, MessageStatus> result = Search.searchWithPeptide(
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
    void searchWithPeptidesFillHitsTest2() throws ParseException {
        List<String> input = new ArrayList<>();
        input.add("IYLGIGLCLLFIVRTLLLHPAIFGLHHIGMQMRIAMFSLIYKKTLKLSSRVLDKISIGQL"); //P13569
        input.add("YVRYFNSSAFFFSGFFVVFLSVLPYALIKGIILRKIFTTISFCIVLRMAVTRQFPWAVQT"); //P13569
        input.add("TDLSQKSLQLESKGLTLNSNAWMNDTVIIDSTVGKDTFFLITWNSLPPSISLWDPSGTIM"); //Q14CN2
        input.add("TGRRGDLATIHGMNRPFLLLMATPLERAQHLQSSRHRRALDTNYCFSSTEKNCCVRQLYI"); //P01137
        input.add("SEWLVLQTPHLEFQEGETIMLRCHSWKDKPLVKVTFFQNGKSQKFSHLDPTFSIPQANHS"); //P12318
        input.add("PPKELVLAGKDAAAEYDELAEPQDFQDDPDIIAFRKANKVGIFIKVTPQREEGEVTVCFK"); //Q9UJW0
        input.add("\tAVLERILAPELSHANATRNLNFSIWNHTPLVLIDERNPHHPMVLDLFGDNHNGLTSSSAS"); //P51168
        input.add("MAPGEKIKAKIKKNLPVTGPQAPTIKELMRWYCLNTNTHGCRRIVVSRGRLRRLLWIGFT"); //P51170
        input.add("FFCNNTTIHGAIRLVCSQHNRMKTAFWAVLWLCTFGMMYWQFGLLFGEYFSYPVSLNINL"); //P37088
        input.add("   AVVENVPPLRWKEFVRRLGLSDHEIDRLELQNGRCLREAQYSMLATWRRRTPRREATLEL"); //P19438
        input.add("ILASPNPDEKTKEELEELMSDIKKTANKVRSKLKSIEQSIEQEEGLNRSSADLRIRKTQH"); //Q16623
        input.add("MTSKLAVALLAAFLISAALCEGAVLPRSAKELRCQCIKTYSKPFHPKFIKELRVIESGPH"); //P10145
        input.add("blabla");

        Pair<List<String[]>, MessageStatus> result = Search.searchWithPeptide(
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
    void searchWithModifiedPeptideInsulinTest() throws IOException, ParseException {

        Pair<List<String[]>, MessageStatus> result = Search.searchWithModifiedPeptide(
                Files.readLines(new File(resourcesPath + "input/ModifiedPeptides/Insulin.txt"), Charset.defaultCharset()),
                MatchType.SUPERSET,
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

        assertEquals(1, hitProteins.size());
        assertEquals(22, hitPathways.size());

        assertEquals(3, iPathways.get("R-HSA-977225").getEntitiesFound().size());
        assertTrue(iPathways.get("R-HSA-977225").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00798:31,00798:43")));
        assertTrue(iPathways.get("R-HSA-977225").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00087:53,00798:31,00798:43")));
        assertFalse(iPathways.get("R-HSA-977225").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308")));
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

    @Test
    void searchWithRsId() {
        HashSet<String> input = new HashSet<>();
        input.add("rs121918101"); //P01308 not found
        input.add("rs28933985"); //P01308 not found
        input.add("rs10840447"); //P01308
        input.add("rs7110099"); //P01308
        input.add("rs555583938"); //P01308

        // Search in the wrong chromosome --> No hits
        imapRsIdsToProteins = (ImmutableSetMultimap<String, String>) getSerializedObject("imapRsIdsToProteins22.gz");
        Search.searchWithRsId(
                input,
                iReactions,
                iPathways,
                imapRsIdsToProteins,
                imapProteinsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                true,
                hitProteins,
                hitPathways
        );

        assertEquals(5, input.size());
        assertEquals(0, hitProteins.size());
        assertEquals(0, hitPathways.size());

        // Search in the right chromosome --> find hits
        imapRsIdsToProteins = (ImmutableSetMultimap<String, String>) getSerializedObject("imapRsIdsToProteins11.gz");
        Search.searchWithRsId(
                input,
                iReactions,
                iPathways,
                imapRsIdsToProteins,
                imapProteinsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                true,
                hitProteins,
                hitPathways
        );

        assertEquals(5, input.size());
        assertEquals(1, hitProteins.size());
        assertTrue(hitProteins.contains("P01308"));
        assertEquals(29, hitPathways.size());
        assertTrue(hitPathways.contains("R-HSA-264876"));
        assertTrue(hitPathways.contains("R-HSA-74749"));
    }

    @Test
    void searchWithChr_Bp() {
        Set<Long> input = new HashSet<>();
        input.add(2176042L); //P01308
        input.add(2176105L); //P01308
        input.add(2176134L); //P01308
        input.add(-1L); //Not found

        // Search in the wrong chromosome --> No hits
        imapChrBpToProteins = (ImmutableSetMultimap<Long, String>) getSerializedObject("imapChrBpToProteins22.gz");
        Search.searchWithChrBp(
                22,
                input,
                iReactions,
                iPathways,
                imapChrBpToProteins,
                imapProteinsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                true,
                hitProteins,
                hitPathways
        );

        assertEquals(4, input.size());
        assertEquals(0, hitProteins.size());
        assertEquals(0, hitPathways.size());

        // Search in the right chromosome --> find hits
        imapChrBpToProteins = (ImmutableSetMultimap<Long, String>) getSerializedObject("imapChrBpToProteins11.gz");
        Search.searchWithChrBp(
                11,
                input,
                iReactions,
                iPathways,
                imapChrBpToProteins,
                imapProteinsToReactions,
                imapReactionsToPathways,
                imapPathwaysToTopLevelPathways,
                true,
                hitProteins,
                hitPathways
        );

        assertEquals(4, input.size());
        assertEquals(1, hitProteins.size());
        assertTrue(hitProteins.contains("P01308"));
        assertEquals(29, hitPathways.size());
        assertTrue(hitPathways.contains("R-HSA-264876"));
        assertTrue(hitPathways.contains("R-HSA-74749"));
    }
}