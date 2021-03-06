package no.uib.pap.methods.search;

import com.google.common.io.Files;
import no.uib.pap.model.InputType;
import no.uib.pap.model.Mapping;
import no.uib.pap.model.MatchType;
import no.uib.pap.model.ProteoformFormat;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.ParseException;

import static org.junit.jupiter.api.Assertions.*;

class SearchModifiedPeptideTest {

    private static String resourcesPath = "../PathwayMatcher/resources/";
    private static String fastaFile = "../PathwayMatcher/resources/uniprot-all.fasta";

    private static Mapping mapping;

    @BeforeAll
    static void loadStaticMapping() {
        mapping = new Mapping(InputType.MODIFIEDPEPTIDE, true);
    }

    @Test
    void singleModifiedPeptideSearchTest() throws IOException, ParseException {
        SearchResult result = Search.searchWithModifiedPeptide(
                Files.readLines(new File(resourcesPath + "input/ModifiedPeptides/SingleModifiedPeptide.txt"), Charset.defaultCharset()),
                mapping,
                true,
                MatchType.SUPERSET,
                0L,
                fastaFile
        );

        assertEquals(1, result.getHitProteins().size());
        assertEquals(11, result.getHitPathways().size());

        assertEquals(2, result.getHitPathwayByStid("R-HSA-168249").getEntitiesFound().size());
        assertTrue(result.getHitPathwayByStid("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;00048:127,00048:132,00048:171,00048:191,00048:226")));
        assertTrue(result.getHitPathwayByStid("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;")));
        assertFalse(result.getHitPathwayByStid("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561;00048:156,00048:161,00048:200,00048:220,00048:255")));
        assertFalse(result.getHitPathwayByStid("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561;00048:200,00048:220")));
    }

    // Warning: When running all tests at once this gets an error, due to shared memory.
    @Test
    void singleModifiedPeptideSearchStrictTest() throws IOException, ParseException {
        SearchResult result = Search.search(
                Files.readLines(new File(resourcesPath + "input/ModifiedPeptides/SingleModifiedPeptide.txt"), Charset.defaultCharset()),
                InputType.MODIFIEDPEPTIDE,
                true,
                mapping,
                MatchType.STRICT,
                0L,
                fastaFile
        );

        assertEquals(1, result.getHitProteins().size());
        assertEquals(11, result.getHitPathways().size());

        assertEquals(1, result.getHitPathwayByStid("R-HSA-168249").getEntitiesFound().size());
        assertTrue(result.getHitPathwayByStid("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;00048:127,00048:132,00048:171,00048:191,00048:226")));
        assertFalse(result.getHitPathwayByStid("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;")));
    }

    @Test
    void singleModifiedPeptideDisplacedSearchTest() throws IOException, ParseException {
        SearchResult result = Search.searchWithModifiedPeptide(
                Files.readLines(new File(resourcesPath + "input/ModifiedPeptides/SingleModifiedPeptideDisplaced.txt"), Charset.defaultCharset()),
                mapping,
                true,
                MatchType.STRICT,
                0L,
                fastaFile
        );

        assertEquals(1, result.getHitProteins().size());
        assertEquals(11, result.getHitPathways().size());

        assertEquals(1, result.getHitPathwayByStid("R-HSA-168249").getEntitiesFound().size());
        assertTrue(result.getHitPathwayByStid("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;00048:127,00048:132,00048:171,00048:191,00048:226")));
        assertFalse(result.getHitPathwayByStid("R-HSA-168249").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("O43561-2;")));
    }

    @Test
    void searchWithModifiedPeptideInsulinTest() throws IOException, ParseException {

        SearchResult result = Search.searchWithModifiedPeptide(
                Files.readLines(new File(resourcesPath + "input/ModifiedPeptides/Insulin.txt"), Charset.defaultCharset()),
                mapping,
                true,
                MatchType.SUPERSET,
                0L,
                fastaFile
        );

        assertEquals(1, result.getHitProteins().size());
        assertEquals(14, result.getHitPathways().size());

        assertEquals(3, result.getHitPathwayByStid("R-HSA-9006934").getEntitiesFound().size());
        assertTrue(result.getHitPathwayByStid("R-HSA-9006934").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00798:31,00798:43")));
        assertTrue(result.getHitPathwayByStid("R-HSA-9006934").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308;00798:95,00798:96,00798:100,00798:109")));
        assertTrue(result.getHitPathwayByStid("R-HSA-9006934").getEntitiesFound().contains(ProteoformFormat.SIMPLE.getProteoform("P01308")));
    }
}