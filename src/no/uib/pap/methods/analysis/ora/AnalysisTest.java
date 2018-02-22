package no.uib.pap.methods.analysis.ora;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSetMultimap;
import no.uib.pap.model.MessageStatus;
import no.uib.pap.model.Pathway;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;

import java.io.InputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

import static org.junit.jupiter.api.Assertions.*;

class AnalysisTest {

    @org.junit.jupiter.api.Test
    void analysis() {
        List<String[]> result = new ArrayList<String[]>();
        MessageStatus messageStatus = new MessageStatus("Sucess", 0, 0, "", "");

        System.out.println(System.getProperty("user.dir"));

        String[] values = {
                "",
                "P01308",
                "R-HSA-1614362",
                "SUMF1 mediates the oxidation of cysteine to formylglycine, producing active arylsulfatases",
                "R-HSA-556833",
                "Metabolism of lipids",
                "R-HSA-1430728",
                "Metabolism"
        };
        result.add(values);

        ImmutableMap<String, Pathway> iPathways = (ImmutableMap<String, Pathway>) getSerializedObject("iPathways.gz");
        ImmutableSetMultimap<String, String> imapProteinsToReactions  = (ImmutableSetMultimap<String,String>) getSerializedObject("imapProteinsToReactions.gz");

        Pair<List<String[]>, MessageStatus> searchResult = new MutablePair<>(result, messageStatus);
        Pair<TreeSet<Pathway>, MessageStatus> analysisResult = Analysis.analysis(iPathways, imapProteinsToReactions.keySet().size(), searchResult.getLeft());

        assertEquals(1, analysisResult.getLeft().size());
        assertTrue(analysisResult.getLeft().contains(iPathways.get("R-HSA-556833")));
        assertEquals(1, analysisResult.getLeft().first().getEntitiesFound().size());
        assertEquals(1, analysisResult.getLeft().first().getReactionsFound().size());
        assertTrue(analysisResult.getLeft().first().getReactionsFound().contains("R-HSA-1614362"));
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