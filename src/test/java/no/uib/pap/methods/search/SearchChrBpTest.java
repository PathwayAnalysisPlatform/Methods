package no.uib.pap.methods.search;

import no.uib.pap.model.InputType;
import no.uib.pap.model.Mapping;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class SearchChrBpTest {

    private static Mapping mapping;

    @BeforeAll
    static void loadStaticMapping() {
        System.out.println("The working directory is: " + System.getProperty("user.dir"));
        mapping = new Mapping(InputType.CHRBP, true);
    }

    @Test
    void searchWithChr_Bp() {
        List<String> input = new ArrayList<>();
        input.add("11 2176042L"); //P01308
        input.add("11 2176105L"); //P01308
        input.add("11 2176134L"); //P01308
        input.add("11 -1L"); //Not found

        SearchResult result = Search.searchWithChrBp(input, mapping, true);

        assertEquals(4, input.size());

        assertEquals(1, result.getHitProteins().size());
        assertTrue(result.getHitProteins().contains("P01308"));
        assertEquals(21, result.getHitPathways().size());
        assertTrue(result.containsPathwayByStid("R-HSA-264876"));
        assertTrue(result.containsPathwayByStid("R-HSA-74749"));
    }
}