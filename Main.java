import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;


public class Main {
    public static void main(String[] args) {


        Graph graph = null;
        try {
            long time;
            //graph = Graph.getInstanceFromEdgeFile(args[0]);

            graph = Graph.getInstanceFromEdgeFile("1000EWG.txt");
            //graph = Graph.getInstance(args);
            //System.out.println(graph.allEdges().size());

            //System.out.println(graph.bridges());


            //System.out.println(graph.allEdges().size());


            time = System.nanoTime();
            //Graph opt = graph.ost();
            System.out.printf("ost took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
            //double optN = opt.numberOfInternalVertices();
            //System.out.println(opt);
            //System.out.println(optN);

            //Graph c = graph.tfpcc();
            //System.out.println(c);


            time = System.nanoTime();
            Graph alg = graph.alg();
            //System.out.println(graph.graphs());

            System.out.printf("alg took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
            double algN = alg.numberOfInternalVertices();


            System.out.println(algN);

            //System.out.println(algN + " / " + optN + " = " + algN / optN);

            PrintWriter pw = null;
            try {
                pw = new PrintWriter("name.txt");
                //頂点を追加するので頂点数分の配列を作ることができず最大のインデックスを代わりに使う
                //pw.println(algN + " / " + optN + " = " + algN / optN);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } finally {
                if (pw != null) pw.close() ;
            }


        } catch (IOException e1) {
            e1.printStackTrace();
        }
    }

}
