import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

public class Graph {

    static IloNumVar[] x;
    static HashMap<String, Integer> map;
    static IloCplex cplex;
    //HashMapよりHashSetの方が良いだろう

    //かならずgetVerticesとgetVertexを使う
    private ArrayList<Vertex> vertices;

    public Graph(double input[][]){
        vertices = new ArrayList<>();
        for(int i = 0;i < input.length;i++){
            Vertex vertex = new Vertex();
            vertex.setIndex(i);
            for(int j = 0; j < input.length;j++){
                if(input[i][j] > 0) {
                    vertex.addEdge(new Edge(i,j));
                }
            }
            vertices.add(vertex);
        }
    }

    public void export(String name){
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(name);
            //頂点を追加するので頂点数分の配列を作ることができず最大のインデックスを代わりに使う
            pw.println(getMaxIndex()+1);
            List<Edge> edges = allEdges();
            pw.println(allEdges().size());
            for(Edge edge: allEdges()){
                pw.format("%d %d %f\n", edge.from, edge.to, 0.0);
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } finally {
            if (pw != null) pw.close() ;
        }
    }

    public Graph(){
        vertices = new ArrayList<>();
    }

    public static Graph getInstanceFromEdgeFile(String fileName) throws IOException {
        Scanner sc = new Scanner(new File(fileName));
        int v = sc.nextInt();
        int e = sc.nextInt();
        double edge[][] = new double[v][v];

        for (int i = 0; i < e; i++) {
            int from = sc.nextInt();
            int to = sc.nextInt();
            double tmp = sc.nextDouble();
            edge[from][to] = 1;
            edge[to][from] = 1;
        }
        return new Graph(edge);
    }

    public static Graph getInstance(String[] args) throws IOException {

         final double COFF = 0.1;

        int num_nodes = 0 ;
        double[][] edge_weight = null;
        Scanner sc = null ;
        PrintWriter pw = null ;

        Graph g;

        try {
            boolean isRandom = true ;
            boolean isEdgeWeighted = true ;
            double sparsity = 0 ;
            if (args.length == 2) {
                sc = new Scanner(System.in) ;
                num_nodes = Integer.parseInt(args[0]) ;
                if (num_nodes < 3) {
                    System.out.println("#vertices must be at least 3.") ;
                    return null ;
                }
                System.out.println("File name: ") ;
                String name = sc.nextLine() ;
                pw = new PrintWriter(name) ;
            } else {
                File f = new File(args[0]) ;
                if (!f.exists()) {
                    System.out.println(args[0] + ": not found.") ;
                    return null ;
                }
                sc = new Scanner(f) ;
                num_nodes = sc.nextInt();
                if (num_nodes < 3) {
                    System.out.println("#vertices must be at least 3.") ;
                    sc.close();
                    return null ;
                }
            }

            //�O���t�̓���
            edge_weight = new double[num_nodes][num_nodes] ;
            if (args.length == 2) {
                Random rand = null ;
                if (isRandom) rand = new Random() ;
                    if (isEdgeWeighted) {
                        for (int i = 0 ; i < num_nodes; i++) {
                            for (int j = i + 1 ; j < num_nodes; j++) {
                                if (!isRandom) {
                                    System.out.println("Weight of edge (" + i + ", " + j + "): ") ;
                                    edge_weight[i][j] = edge_weight[j][i] = sc.nextInt();
                                }else {
                                    setEdgeWeight(num_nodes, edge_weight, COFF);
                                }
                            }
                        }
                    } else if (!isRandom){
                        System.out.println("Input the edges! Quit by inputting a negative endpoint!");
                        while (true) {
                            System.out.println("From: ") ;
                            int from = sc.nextInt() ;
                            System.out.println("To: ") ;
                            int to = sc.nextInt() ;
                            if (from < 0 || to < 0) break ;
                            edge_weight[from][to] = edge_weight[to][from] = 1;
                        }
                    } else {
                        for (int i = 0 ; i < num_nodes; i++) {
                            for (int j = i + 1 ; j < num_nodes; j++) {
                                if ( rand.nextDouble() <= sparsity ) {
                                    edge_weight[i][j] = edge_weight[j][i] = 1;
                                }
                            }
                        }
                    }

                for (int i = 0 ; i < num_nodes; i++)
                    edge_weight[i][i] = 0 ;


            } else {
                for (int i = 0 ; i < num_nodes; i++) {
                    for (int j = 0 ; j < num_nodes; j++)
                        edge_weight[i][j] = sc.nextDouble( ) ;
                }

            }

            g = new Graph(edge_weight);
            /*
            while(g.graphs().size() > 1){
                setEdgeWeight(num_nodes,edge_weight,COFF);
                g = new Graph(edge_weight);
            }
            */

            if (pw != null) {//save the instance to file
                pw.println(num_nodes) ;
                for (int i = 0 ; i < num_nodes; i++) {
                    for (int j = 0 ; j < num_nodes; j++)
                        pw.print(edge_weight[i][j] + " ") ;
                    pw.println("") ;
                }

            }
        } finally {
            if (sc != null) sc.close() ;
            if (pw != null) pw.close() ;
        }

        //�O���t������ĕԂ��D

        System.out.println("instance created");
        return g;
    }

    static void setEdgeWeight(int num_nodes, double[][] edge_weight, double coff){
        Random rand = new Random() ;
        for (int i = 0 ; i < num_nodes; i++) {
            for (int j = i + 1 ; j < num_nodes; j++) {
                    edge_weight[i][j] = edge_weight[j][i] =
                            rand.nextDouble() < coff ? 1 : 0;
            }
        }
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for(Vertex vertex: vertices) {
            sb.append("Vertex: " + vertex.getIndex() + "\n");
            for(Edge edge: vertex.getEdges()){
                sb.append("Edge: " + edge.toString() + ", \n");
            }
        }
        return sb.toString();
    }

    public ArrayList<Vertex> getVertices() {
        return vertices;
    }

    public Vertex getVertex(int index) {
        for(Vertex vertex: vertices) {
            if(vertex.getIndex() == index) return vertex;
        }
        return null;
    }

    //edge.toからverticesのindexがほしいときに使う
    //TODO: 名前紛らわしい
    public int getVertexIndex(int index) {

        for(int i = 0; i < vertices.size(); i++){
            if(vertices.get(i).getIndex() == index) return i;
        }

        return -1;
    }

    public Vertex deleteVertex(int index) {

        for(Iterator<Vertex> it = vertices.iterator(); it.hasNext(); ) {
            Vertex vertex = it.next();
            if(vertex.getIndex()== index) {

                //帰る方の辺が削除されてないため
                for(Edge edge : vertex.getEdges()){
                    deleteOnlySwappedEdge(edge);
                }

                it.remove();
                return vertex;
            }
        }
        return null;
    }

    public void deleteVertex(Vertex vertex) {
        //帰る方の辺が削除されてないため
        for(Edge edge : vertex.getEdges()){
            deleteOnlySwappedEdge(edge);
        }
        vertices.remove(vertex);
    }



    public void deleteEdge(Edge edge) {
        edges(edge.from).removeIf(edge1 -> edge1.to == edge.to);
        edges(edge.to).removeIf(edge1 -> edge1.to == edge.from);
    }

    //頂点削除の際にdeleteEdgeを実行するとiterator実行中のエラーになるため
    private void deleteOnlySwappedEdge(Edge edge){
        //vertices.get(edge.to).getEdges().remove(edge.swapEdge())
        Vertex v = getVertex(edge.to);
        if(v != null) edges(edge.to).removeIf(edge1 -> edge1.to == edge.from);
    }

    public LinkedList<Boolean> cycleChecks(){
        LinkedList<Boolean> ret = new LinkedList<>();
        for (int i = 0; i < vertices.size(); i++) {
            boolean isVisited[] = new boolean[vertices.size()];
            ret.add(cycleCheck(i,-1,isVisited));
        }
        return ret;
    }

    public boolean isCycle(){
        //System.out.println("");
        boolean isVisited[] = new boolean[vertices.size()];
        return cycleCheck(0,-2,isVisited);

    }

    boolean cycleCheck(int start,int prev, boolean[] isVisited){
        isVisited[start] = true;
        LinkedList<Edge> edges = vertices.get(start).getEdges();
        //if(edges.size() > 2) return false;
        for(Edge edge: edges){
            int toIndex = getVertexIndex(edge.to);
            if(toIndex == prev) continue;
            if(toIndex == -1) continue;//TODO
            //System.out.println("to: " + edge.to +  ", index: " + toIndex);
            if(isVisited[toIndex] || cycleCheck(toIndex, start, isVisited)){
                return true;
            }
        }
        return false;
    }

    public Graph(int n){
        vertices = new ArrayList<Vertex>(n);
        for(int i = 0; i < n; i++){
            vertices.add(new Vertex(i));
        }
    }

    public Graph(LinkedList<Integer> component) {
        vertices = new ArrayList<Vertex>(component.size());
        for(Integer i: component){
            vertices.add(new Vertex(i));
        }
    }

    private int getMaxIndex(){
        int max = -1;
        for(Vertex vertex: vertices){
            if(vertex.getIndex() > max) max = vertex.getIndex();
        }
        return max;
    }

    //新しく頂点を追加し，その番号を返す
    //頂点と辺を追加したい場合は
    //危なっかしいのでとりあえず使うのやめよう
    /*
    public int addVertex(){
        int nextIndex = getMaxIndex() + 1;
        vertices.add(new Vertex(nextIndex));
        return nextIndex;
    }
    */

    public int addVertex(Vertex vertex){
/*
        for(Vertex v: vertices){
            if(v.getIndex() == vertex.getIndex()) return vertex.getIndex();
        }
        */
        if(!vertices.contains(vertex)) vertices.add(vertex);
        for(Edge edge: vertex.getEdges()) {
            if(getVertex(edge.to) == null) addVertex(new Vertex(edge.to));
            addOnlySwappedEdge(edge); //頂点を戻した際に周りの辺から行く辺を追加する必要がある
        }
        return vertex.getIndex(); //必要ないかも
    }

    public LinkedList<Edge> edges(int i){
        //vertexがないときnullが変えるので注意
        return getVertex(i).getEdges();
    }



    public void addEdge(Edge edge) {

        Vertex fromV = getVertex(edge.from);
        if(fromV == null) {
            fromV = new Vertex(edge.from);
            addVertex(fromV);
        }

        for(Edge e: fromV.getEdges()){
            if(e.to == edge.to){
                //System.out.println("同じ辺を追加しようとしています");
                return;
            }
        }

        fromV.addEdge(edge);

/*
        if(!fromV.containsEdge(edge.to)){
            fromV.addEdge(edge);
        }
        */

        Vertex toV = getVertex(edge.to);
        if(toV == null) {
            toV = new Vertex(edge.to);
            addVertex(toV);
        }

        for(Edge e: toV.getEdges()){
            if(e.to == edge.from) {
                //System.out.println("同じ辺を追加しようとしています2");
                return;
            }
        }

        toV.addEdge(edge.swapEdge());
/*
        Edge rEdge = edge.swapEdge();
        if(!toV.containsEdge(rEdge.to)){
            toV.addEdge(rEdge);
        }
*/
    }

    public void addOnlySwappedEdge(Edge edge) {
        edges(edge.to).add(edge.swapEdge());
    }


    public List<Edge> allEdges(){
        long time = System.nanoTime();
        List<Edge> ret = new LinkedList<>();
        for(Vertex vertex: vertices){
            for(Edge edge: vertex.getEdges()){
                if(edge.from < edge.to){
                    ret.add(new Edge(edge.from, edge.to));
                }
            }
        }

        //System.out.printf("allEdge took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
        return ret;
    }

    private Graph fromCplex(IloNumVar[] x, double[] val){
        Graph newGraph = new Graph();

        for (int j = 0; j < val.length; ++j){
            if(val[j] > 0) {
                Matcher m =  Pattern.compile("_([0-9]+)_([0-9]+)").matcher(x[j].getName());
                if(m.find()){
                    int from = Integer.parseInt(m.group(1));
                    int to = Integer.parseInt(m.group(2));
                    if(newGraph.getVertex(from) == null) newGraph.addVertex(new Vertex(from));
                    if(newGraph.getVertex(to) == null) newGraph.addVertex(new Vertex(to));
                    newGraph.addEdge(new Edge(from, to));
                }
            }
        }
        return newGraph;
    }

    private boolean existsEdge(int from, int to){
        if(getVertex(from) == null) return false;
        if(getVertex(to) == null) return false;
        return getVertex(from).containsEdge(to);
    }

    public LinkedList<LinkedList<Integer>> getConnectedComponents() {


        int nodeNum = vertices.size();

        boolean isVisited[] = new boolean[nodeNum];

        LinkedList<LinkedList<Integer>> components = new LinkedList<>();
        boolean addflag = false;
        loop1: for(int i = 0; i < nodeNum; i++){

            for(LinkedList set: components){
                if(set.contains(i) && !components.isEmpty()) continue loop1;
            }
            LinkedList<Integer> component = new LinkedList<>();
            LinkedList<Integer> q = new LinkedList<>();
            component.add(i);
            q.add(i);
            isVisited[i] = true;
            while(!q.isEmpty()){
                int w = q.poll();
                for(Edge edge: edges(w)){
                    if(isVisited[edge.to]) continue;

                    isVisited[edge.to] = true;
                    component.add(edge.to);
                    q.add(edge.to);

                }
            }
            components.add(component);
        }
        return components;
    }

    //componentに乗っかる形なのでできればgraphのみで判断をしたい
    public LinkedList<Graph> graphs() {
        int nodeNum = vertices.size();
        boolean isVisited[] = new boolean[nodeNum];
        LinkedList<Graph> graphs = new LinkedList<>();

        LinkedList<LinkedList<Integer>> components = new LinkedList<>();
        loop1: for(int i = 0; i < nodeNum; i++){

            for(LinkedList component: components){
                if(component.contains(i) && !components.isEmpty()) continue loop1;
            }
            LinkedList<Integer> component = new LinkedList<>();// indexを覚える
            LinkedList<Integer> q = new LinkedList<>(); //indexを覚える
            Graph graph = new Graph();
            component.add(i);
            q.add(i);
            graph.vertices.add(vertices.get(i));
            isVisited[i] = true;
            while(!q.isEmpty()){
                int fromIndex = q.poll();
                Vertex vertex = vertices.get(fromIndex);
                for(Edge edge: vertex.getEdges()){
                    int toIndex = getVertexIndex(edge.to);
                    if(isVisited[toIndex]) continue;

                    isVisited[toIndex] = true;
                    component.add(toIndex);
                    q.add(toIndex);
                    //0を追加した際に辺の先の頂点を既に追加済みになっているためaddVertexを使うと重複追加になる
                    graph.vertices.add(vertices.get(toIndex));
                }
            }
            //System.out.println(graph);
            graphs.add(graph);
            components.add(component);
        }
        return graphs;
    }

    public LinkedList<Graph> graphs2() {
        int nodeNum = vertices.size();
        boolean isVisited[] = new boolean[nodeNum];
        LinkedList<Graph> graphs = new LinkedList<>();
        Graph graph;
        for (int i = 0; i < nodeNum; i++) {
            if(!isVisited[i]) {
                isVisited[i] = true;
                graph = new Graph();
                graph.addVertex(new Vertex(vertices.get(i).getIndex()));
                graphsDFS(i, isVisited, graph);
                for(Edge edge: vertices.get(i).getEdges()){
                    graph.addEdge(new Edge(edge.from, edge.to));
                }
                graphs.add(graph);
            }
        }
        return graphs;
    }

    private void graphsDFS(int current, boolean[] isVisited, Graph graph){
        Vertex v = vertices.get(current);
        for(Edge edge: v.getEdges()){
            int toI = getVertexIndex(edge.to);
            Edge newEdge = new Edge(current, edge.to);
            if(isVisited[toI]){
                //if(!v.getEdges().contains(newEdge)) graph.addEdge(newEdge);
                continue;
            }
            isVisited[toI] = true;
            graph.addVertex(new Vertex(v.getIndex()));
            //graph.addVertex(new Vertex(edge.to));
            graphsDFS(toI, isVisited, graph);
            //graph.addEdge(newEdge);
        }
    }



    public LinkedList<Edge> bridges(){
        LinkedList<Edge> bridges = new LinkedList<>();


        //edgeのtoをそのまま配列の添字に使わないこと
        int size = vertices.size();

        if(size == 0) return new LinkedList<>();

        int low[] = new int[size];
        Arrays.fill(low, -1);

        int pre[] = new int[size];
        Arrays.fill(pre, -1);

        int count[] = {0};

        bridgesDFS(bridges, 0, count, low, pre);

        return bridges;
    }

    private int bridgesDFS(LinkedList<Edge> bridges, int from , int[] count , int low[] ,int pre[]) {
        count[0]++;
        pre[from] = count[0];
        low[from] = pre[from];

        //配列の添字に合わせるためgetvを使う
        for(Edge edge: vertices.get(from).getEdges()){

            int to = getVertexIndex(edge.to);

            //未訪問の場合
            if(pre[to] == -1) {

                //fromの最小値を子の最小値と比較し採用する
                low[from] = Math.min(low[from], bridgesDFS(bridges, to, count, low, pre));

                //もし他にfromとtoを繋ぐルートがなければ追加a
                if(low[to] == pre[to]) {
                    //from < toとなるように調整し，橋に追加
                    if(edge.from < edge.to) bridges.add(edge);
                    else bridges.add(new Edge(edge.to, edge.from));
                }
            } else {
                //訪問済みであれば最小値を更新
                low[from] = Math.min(low[from], low[to]);
            }
        }
        return low[from];
    }

    public Graph opr1() {
        if(vertices.size() <= 3) return this;
        long time = System.nanoTime();
        LinkedList<Vertex> leaves = new LinkedList<Vertex>();
        int lSize = 0;

        for(Vertex vertex: getVertices()) {
            if(vertex.getEdges().size() == 1) {
                leaves.add(vertex);
                lSize++;
            }
        }

        for (int i = 0; i < lSize - 1; i++) {
            for (int j = i + 1; j < lSize; j++) {
                Vertex v1 = leaves.get(i);
                Vertex v2 = leaves.get(j);
                if(v1.getEdges().getFirst().to == v2.getEdges().getFirst().to) {
                    System.out.println("opr1 deleted:" + v1);
                    deleteVertex(v1);
                    break;
                }
            }
        }


        System.out.printf("opr1 took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
        return this;
    }


    public Graph opr2() {
        System.out.println("---OPR2 START---");
        long time = System.nanoTime();


        int vSize= vertices.size();

        //u1もu2もカットポイント
        List<Edge> allEdges = allEdges();
        TreeSet<Integer> cutPoints = cutPoints();
        int fromI= -1;
        int fromSize = -1;
        int toI= -1;
        int toSize = -1;
        for(Edge edge: allEdges) {

            //キャッシュ
            if(edge.from != fromI){
                fromI = edge.from;
                fromSize = getVertex(edge.from).getEdges().size();
            }

            if(edge.to != toI){
                toI = edge.to;
                toSize = getVertex(edge.to).getEdges().size();
            }


            //LinkedList<Edge> bridges  = bridges();
            //連結成分Kを作るには3つは辺が必要
            if(fromSize < 3 || getVertex(edge.to).getEdges().size() < 3){
            //if(bridges.contains(edge)) {
                System.out.println("continue");
                continue;
            }

            //System.out.println(edge);
            if(!cutPoints.contains(edge.from) || !cutPoints.contains(edge.to)) continue;

            if(fromSize < 3 || toSize < 3) continue;
            if(opr2CheckComponent(edge, vSize) && opr2CheckComponent(edge.swapEdge(),vSize)) {
                System.out.println("opr2 deleted: " + edge);
                deleteEdge(edge);
            }

        }

        System.out.printf("opr2 took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
        return this;
    }

    //ボトルネック
    private boolean opr2CheckComponent(Edge edge, int vSize) {
        long time = System.nanoTime();
        List<Edge> uEdges = new LinkedList<>();
       // int maxPlus1 = getMaxIndex() + 1;
        for(Edge e: edges(edge.from)){
            uEdges.add(new Edge(e.from, e.to));
        }

        //uEdges = edges(edge.from);

        Vertex deletedVertex = deleteVertex(edge.from);


        boolean[] flag = new boolean[vSize -1];
        assinCol(getVertexIndex(edge.to),flag);
        for(Edge u1Edge: uEdges){
            boolean isSame = flag[getVertexIndex(u1Edge.to)];
            if(!isSame){
                addVertex(deletedVertex);
                System.out.printf("opr2check took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
                return true;
            }
        }


        addVertex(deletedVertex);
        System.out.printf("opr2check took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
        return false;
    }

    //TODO: flagを返してあげて何回でもO(1)で確認できるようにする
    public void isSameComponent(int fromI, int to, boolean flag[]){

        int toI = getVertexIndex(to);
        flag[fromI] = true;
        assignDFS(fromI, flag);

    }

    public void assinCol(int fromI, boolean flag[]){
        flag[fromI] = true;
        assignDFS(fromI, flag);
    }

    public void assignDFS(int current, boolean[] flag){
        boolean ret = false;
        for(Edge edge: vertices.get(current).getEdges()){
            int toI = getVertexIndex(edge.to);
            if(flag[toI]) continue;
            flag[toI] = true;
            assignDFS(toI, flag);
        }
    }

    //
    public UnionFind createUnionFind(){
        //辺一覧を取ってくるAllEdgesが欲しい
        int n = vertices.size();
        UnionFind uf = new UnionFind(n);
        for(Edge edge: allEdges()){
            uf.unite(getVertexIndex(edge.from), getVertexIndex(edge.to));
        }
        return uf;
    }

    public UnionFind createUnionFind(TreeSet<Edge> edges){
        int n = vertices.size();
        UnionFind uf = new UnionFind(n);
        for(Edge edge: edges){
            uf.unite(getVertexIndex(edge.from), getVertexIndex(edge.to));
        }
        return uf;
    }


    public Graph opr3() {
        System.out.println("---OPR3 START---");
        long time = System.nanoTime();

        for(Edge bridge: bridges()){
            //もしブリッジの端が１つの点だった場合バグる
            if(getVertex(bridge.to).getEdges().size() == 1 || getVertex(bridge.from).getEdges().size() == 1) continue;
            deleteEdge(bridge);

            //ブリッジの場合，連結成分が増えるか
            LinkedList<Graph> graphs = graphs();
            Graph g1 = graphs.get(0);
            Graph g2 = graphs.get(1);
            Vertex u1, u2;
            if(g1.getVertexIndex(bridge.from) >= 0) {
                u1 = getVertex(bridge.from);
                u2 = getVertex(bridge.to);
            } else {
                u1 = getVertex(bridge.to);
                u2 = getVertex(bridge.from);
            }


            if(g1.cutPoints().contains(u1.getIndex()) && g2.cutPoints().contains(u2.getIndex())) {

                System.out.printf("opr3 took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
                Graph merged = g1.alg().merge(g2.alg());
                merged.addEdge(bridge);
                return merged;
            } else {
                addEdge(bridge);
            }

        }
        //なにもなければそのまま
        System.out.println("opr3 nothing");
        return this;
    }

    public TreeSet<Integer> cutPoints() {

        TreeSet<Integer> cutPoints = new TreeSet(); //求めた切断点

        int size = vertices.size();

        boolean visited[] = new boolean[size];
        int pre[] = new int[size]; //頂点を発見した順番を覚える
        int low[] = new int[size]; //各頂点における最小値
        int parent[] = new int[size]; //その頂点の親になる点
        boolean ap[] = new boolean[size]; // その頂点が切断点であるかどうか

        int[] count = {0}; //何番目の訪問かカウント

        // 初期化
        for (int i = 0; i < size; i++) {
            parent[i] = -1;
            visited[i] = false;
            ap[i] = false;
        }

        // 各頂点を始点としてDFSを実行(連結じゃない場合でも動くようにするため)
        for (int i = 0; i < size; i++)
            if (!visited[i])
                cutPointsDFS(i, visited, pre, low, parent, ap, count);

        // 切断点を戻り値に含める
        for (int i = 0; i < size; i++)
            if (ap[i])
                cutPoints.add(vertices.get(i).getIndex());

        return cutPoints;
    }

    void cutPointsDFS(int from, boolean visited[], int pre[], int low[], int parent[], boolean ap[], int[] count) {

        // 子の数を覚える
        int children = 0;

        //頂点を訪問済に
        visited[from] = true;

        // 訪問順および最小値を初期化
        count[0]++;
        pre[from] = low[from] = count[0];


        for(Edge edge: vertices.get(from).getEdges()) {

            //行き先のインデックスを取得(配列に合わせる)
            int to = getVertexIndex(edge.to);

            // もしtoが訪問済みでなければ子に追加し，fromを親にセットして再帰
            if (!visited[to]) {
                children++;
                parent[to] = from;
                cutPointsDFS(to, visited, pre, low, parent, ap, count);

                // もしfromの最小値の方が小さければfromも更新
                low[from]  = Math.min(low[from], low[to]);


                // fromが根であり，子が2つ以上である場合，切断点
                if (parent[from] == -1 && children > 1)
                    ap[from] = true;

                // fromが根ではないが，fromの子からtoに向かう辺がない場合，切断点
                if (parent[from] != -1 && low[to] >= pre[from])
                    ap[from] = true;
            }

            // 最小値を更新
            else if (to != parent[from])
                low[from]  = Math.min(low[from], pre[to]);
        }
    }

    public Graph opr4(){

        long time = System.nanoTime();
        System.out.println("opr4 start");
        //cut pointが最後の点の場合，先に削除されているため正しく+1されていない問題が起きたため外に出す
        int nextIndex = getMaxIndex() + 1;
        for(Integer i: cutPoints()){

            //cutPointsは中身のindexが反映されているため
            Vertex deletedVertex = deleteVertex(i);

            System.out.println("opr4 cutPoint: " + i);
            boolean flag = false;
            for(Graph graph: graphs()){
                int size = graph.vertices.size();
                LinkedList<Vertex> deletes = new LinkedList<>();
                if(2 <= size && size <= 8) {
                    LinkedList<Integer> deletedComponent = new LinkedList<>();
                    //グラフから連結成分の点を削除し，頂点を戻し，点を追加しその辺を追加す
                    //再生成してない
                    Graph k = graph.alg();
                    for (Vertex vertex : graph.vertices) {
                        System.out.println("opr4 deleted:" + vertex);
                        deletedComponent.add(vertex.getIndex());
                        deleteVertex(vertex);
                    }

                    System.out.println("opr4 readd: " + deletedVertex);
                    Edge newEdge = new Edge(deletedVertex.getIndex(), nextIndex);

                    Vertex v = new Vertex(deletedVertex.getIndex());

                    //削除する成分以外の辺を追加
                    for (Edge edge : deletedVertex.getEdges()) {
                        if (deletedComponent.contains(edge.to)) continue;
                        addEdge(new Edge(v.getIndex(), edge.to));
                    }


                    Vertex newVertex = new Vertex(nextIndex);

                    System.out.println("opr4 new added: " + newVertex);
                    addVertex(newVertex);
                    System.out.println("opr4 add edge: " + newEdge);
                    addEdge(newEdge);


                    System.out.printf("opr4 took %.3f ms to run%n", (System.nanoTime() - time) / 1e6);
                    Graph g1 = alg();
                    g1.deleteVertex(newVertex);
                    g1 = g1.merge(k);

                    for (Edge edge : deletedVertex.getEdges()) {
                        if (deletedComponent.contains(edge.to)) {
                            g1.addEdge(new Edge(v.getIndex(), edge.to));
                        }

                    }
                    return g1;
                }
            }

            //操作を行えない場合は元に戻す
            addVertex(deletedVertex);
        }

        //操作を行えない場合のみ呼ばれる
        return this;
    }

    public Graph merge(Graph graph) {
        Graph newGraph = new Graph();
        for(Vertex vertex: vertices){
            newGraph.addVertex(new Vertex(vertex.getIndex()));
        }

        for(Vertex vertex: vertices){
            for(Edge edge: vertex.getEdges()){
                newGraph.addEdge(new Edge(edge.from, edge.to));
            }
        }

        for(Vertex vertex: graph.vertices){
            newGraph.addVertex(new Vertex(vertex.getIndex()));
        }

        for(Vertex vertex: graph.vertices){
            for(Edge edge: vertex.getEdges()){
                newGraph.addEdge(new Edge(edge.from, edge.to));
            }
        }

        return newGraph;
    }

    public Graph alg(){
        Graph g = opr1().opr2().opr3().opr4();
        if(g.vertices.size() <= 8) return g.ost();
        Graph c = g.tfpcc();


        Graph preprocessedC = c.opr5(g).opr6(g).opr7(g);
        return preprocessedC.finalOpr1(g).finalOpr2(g).finalOpr3(g);
        //return c;
    }


    //三角形を除く条件を後から追加する必要がある
    public Graph tfpcc(){
        long time = System.nanoTime();
        try {
            cplex = new IloCplex();
            cplex.setParam(IloCplex.IntParam.Threads, Runtime.getRuntime().availableProcessors() - 2);
            map = new HashMap<>();

            int edgesSize = allEdges().size();

            x = cplex.boolVarArray(edgesSize);
            int count = 0;
            for(Edge edge: allEdges()){
                //変数名と同じで最初は数字を入れられない
                String key = "_" + edge.from + "_" + edge.to;
                map.put(key, count);
                x[count].setName(key);
                count++;
            }

            //各頂点からは1本か2本
            for(Vertex vertex: vertices){
                IloLinearNumExpr constraint = cplex.linearNumExpr();
                for(Edge edge: vertex.getEdges()){
                    String key;

                    if(edge.from < edge.to) key = toKey(edge.from, edge.to);
                    else                    key = toKey(edge.to, edge.from);
                    constraint.addTerm(1, x[map.get(key)]);
                }

                cplex.addLe(1, constraint);
                cplex.addGe(2, constraint);
            }

            //最大化
            IloLinearNumExpr maximizeFunc = cplex.linearNumExpr();
            for(Edge edge: allEdges()){
                maximizeFunc.addTerm(1, x[map.get(toKey(edge.from, edge.to))]);
            }
            cplex.addMaximize(maximizeFunc);



            cplex.use(new TriangleCutCallback());



            cplex.exportModel("tfpcc.lp");
            if (cplex.solve()) {
                //cplex.output().println("Solution status = " + cplex.getStatus());
                //cplex.output().println("Solution value = " + cplex.getObjValue());

                Graph ret = fromCplex(x, cplex.getValues(x));
                //System.out.println("RET::: " + ret);
                System.out.printf("TFPCC took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
                return ret;

            }
            cplex.end();

        } catch (IloException e) {
            e.printStackTrace();
        }
        return null;
    }

    static private class TriangleCutCallback extends IloCplex.LazyConstraintCallback {

        public TriangleCutCallback() {}

        @Override
        protected void main() throws IloException {
            Graph tmpGraph = new Graph();
            Graph g = tmpGraph.fromCplex(x, this.getValues(x));

            LinkedList<Graph> triangles = new LinkedList<>();

            for(Graph graph: g.graphs()){
                if(graph.vertices.size() == 3 && graph.isCycle()) {
                    triangles.add(graph);
                }
            }

            for(Graph triagle : triangles){
                IloLinearNumExpr over1 = cplex.linearNumExpr();
                for(Vertex tv: triagle.getVertices()) {
                    for(Vertex v: g.getVertices()){

                    //元のGに辺のあるもののみだけ追加
                        String key;
                        int from, to;
                        if(tv.getIndex() < v.getIndex()){
                            from = tv.getIndex();
                            to = v.getIndex();
                        } else {
                            from = v.getIndex();
                            to = tv.getIndex();
                        }

                        key = _toKey(from, to);
                        //System.out.println(map);
                        if(map.containsKey(key)) {
                            //三角形の辺を削除
                            if(triagle.existsEdge(from, to)){
                                continue;
                            }
                            //そこ以外の辺で条件追加
                            //System.out.println(key);
                            over1.addTerm(1, x[map.get(key)]);
                        }
                    }

                }
                //ここで1以上の条件追加
                this.add(cplex.ge(over1,1));
            }
        }
    }

    @Override
    protected Graph clone() {
        Graph copied = new Graph();
        for(Vertex vertex: vertices) {
            Vertex newVertex = new Vertex(vertex.getIndex());

            for(Edge edge: vertex.getEdges()){
                newVertex.addEdge(new Edge(edge.from, edge.to));
            }

            copied.addVertex(newVertex);
        }
        return copied;
    }

    public Graph ost(){
        long time = System.nanoTime();
        try {
            cplex = new IloCplex();
            map = new HashMap<>();
            //IloCplex.Param.Threads(1);
            cplex.setParam(IloCplex.IntParam.Threads, Runtime.getRuntime().availableProcessors() - 2);

            int edgesSize = allEdges().size();
            int verticesSize = vertices.size();

            x = cplex.boolVarArray(edgesSize);
            int count = 0;
            for(Edge edge: allEdges()){
                //変数名と同じで最初は数字を入れられない
                String key = toKey(edge.from, edge.to);
                map.put(key, count);
                x[count].setName(key);
                count++;
            }

            //辺の数は頂点数-1
            IloLinearNumExpr edgesSizeConstraint = cplex.linearNumExpr();
            for(Edge edge : allEdges()){
                edgesSizeConstraint.addTerm(1.0, x[map.get(toKey(edge.from, edge.to))]);
            }
            cplex.addEq(verticesSize - 1, edgesSizeConstraint);

/*
            //各頂点から出ている辺は1本以上
            for(Vertex vertex: vertices){
                IloLinearNumExpr constraint = cplex.linearNumExpr();
                for(Edge edge: vertex.getEdges()){
                    String key;
                    if(edge.from < edge.to){
                        key = toKey(edge.from, edge.to);
                    } else {
                        key = toKey(edge.to, edge.from);
                    }
                    constraint.addTerm(1, x[map.get(key)]);

                }
                cplex.addLe(1, constraint);
            }
*/



            IloNumVar[] v = cplex.boolVarArray(verticesSize);
            //IloNumVar[] v = cplex.numVarArray(verticesSize,0,1);
            int vi = 0;
            for(Vertex vertex: vertices){
                v[vi].setName("_" + vertex.getIndex());
                IloLinearNumExpr constraint = cplex.linearNumExpr();
                for(Edge edge: vertex.getEdges()){
                    String key;
                    if(edge.from < edge.to){
                        key = toKey(edge.from, edge.to);
                    } else {
                        key = toKey(edge.to, edge.from);
                    }
                    constraint.addTerm(1, x[map.get(key)]);

                }
                //constraint.addTerm(1 - vertex.getEdges().size(), v[vi]);
                constraint.addTerm(-1, v[vi]);
                cplex.addLe(1, constraint);
                vi++;
            }


            IloLinearNumExpr maximizeFunc = cplex.linearNumExpr();
            for(int i = 0; i < v.length;i++){
                maximizeFunc.addTerm(1, v[i]);
            }

            cplex.addMaximize(maximizeFunc);

            //連結成分は1つという条件が必要
            cplex.use(new CycleCutCallback());
            cplex.exportModel("ost.lp");

            if (cplex.solve()) {
                //cplex.output().println("Solution status = " + cplex.getStatus());
                //cplex.output().println("Solution value = " + cplex.getObjValue());

                Graph ret = fromCplex(x, cplex.getValues(x));
                System.out.println("OST: " + ret);
                System.out.printf("OST took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
                return ret;


            }
            cplex.end();

        } catch (IloException e) {
            e.printStackTrace();
        }

        return null;

    }


    static private class CycleCutCallback extends IloCplex.LazyConstraintCallback {

        public CycleCutCallback() {}

        @Override
        protected void main() throws IloException {
            System.out.println("aoao");
            Graph tmpGraph = new Graph();
            Graph g = tmpGraph.fromCplex(x, this.getValues(x));

            LinkedList<Graph> components = new LinkedList<>();

            LinkedList<Graph> graphs = g.graphs();

            if(graphs.size() > 1){

                components.addAll(graphs);

                for(Graph component : components){
                    IloLinearNumExpr over1 = cplex.linearNumExpr();
                    for(Vertex tv: component.getVertices()) {
                        for(Vertex v: g.getVertices()){

                            //元のGに辺のあるもののみだけ追加(mapを使う)
                            String key;
                            int from, to;
                            if(tv.getIndex() < v.getIndex()){
                                from = tv.getIndex();
                                to = v.getIndex();
                            } else {
                                from = v.getIndex();
                                to = tv.getIndex();
                            }

                            key = _toKey(from, to);
                            //System.out.println(map);
                            if(map.containsKey(key)) {
                                //三角形の辺を削除
                                if(component.existsEdge(from, to)){
                                    continue;
                                }
                                //そこ以外の辺で条件追加
                                //System.out.println(key);
                                over1.addTerm(1, x[map.get(key)]);
                            }
                        }

                    }
                    //ここで1以上の条件追加
                    this.add(cplex.ge(over1,1));
                }
            }
        }
    }

    private Graph finalOpr1(Graph g){
        long time = System.nanoTime();
        LinkedList<Graph> graphs =  graphs();
        ArrayList<Graph> cycles = new ArrayList<>();
        for(Graph graph: graphs) {
            if(graph.isCycle()) {
                cycles.add(graph);
            }
        }

        int n = cycles.size();

        if(n < 2) {


            return this;
        }

        label: for(int i = 0;i < n - 1;i++){
            Graph c1 = cycles.get(i);


            for(int j = i + 1;j < n;j++){

                Graph c2 = cycles.get(j);

                for(Vertex v1: c1.vertices){

                    //c1の頂点をGでみたとき，辺を見ていったときにc2の頂点が見つかるか
                    for(Edge edge: g.getVertex(v1.getIndex()).getEdges()) {
                        Vertex v2 = c2.getVertex(edge.to);
                        if(v2 != null) {
                            Edge deleted1 = getVertex(v1.getIndex()).getEdges().getFirst();
                            Edge deleted2 = getVertex(v2.getIndex()).getEdges().getFirst();
                            deleteEdge(deleted1);
                            deleteEdge(deleted2);
                            addEdge(edge);
                            continue label;
                        }
                    }

                }
            }
        }

        System.out.printf("final1 took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
        return this;
    }
    public Graph finalOpr2(Graph g){
        long time = System.nanoTime();
        UnionFind uf = createUnionFind();
        LinkedList<Graph> graphs = graphs();

        //export("c.txt");
        //g.export("g.txt");


        for(Graph graph: graphs) {

            //接続先が弄られるのに注意
            if(graph.isCycle()) {
                label: for(Vertex vertex: graph.vertices) {
                    //Gにおいての頂点を取ってくる
                    Vertex vertexOnG =  g.getVertex(vertex.getIndex());
                    for(Edge edgeOnG: vertexOnG.getEdges()) {
                        //気をつける
                        int fromIndex = getVertexIndex(edgeOnG.from);
                        int toIndex = getVertexIndex(edgeOnG.to);

                        if(!uf.same(fromIndex, toIndex)) {
                            Edge deleted = vertex.getEdges().getFirst();
                            deleteEdge(deleted);
                            addEdge(edgeOnG);
                            uf.unite(fromIndex, toIndex);
                            break label; //一度実行するとサイクルではなくなるため
                        }
                    }
                }
            }
        }

        //graphsには問題なし，しかしisCycleでバグるということは以下のどこかで書き換えが発生
        /*
        for(Graph graph: graphs) {

            if(graph.isCycle()) {
                for(Vertex vertex: vertices) {
                    //Gにおいての頂点を取ってくる
                    Vertex vertexOnG =  g.getVertex(vertex.getIndex());
                    for(Edge edgeOnG: vertexOnG.getEdges()) {
                        //気をつける
                        int fromIndex = getVertexIndex(edgeOnG.from);
                        int toIndex = getVertexIndex(edgeOnG.to);

                        if(!uf.same(fromIndex, toIndex)) {
                            Edge deleted = vertex.getEdges().getFirst();
                            deleteEdge(deleted);
                            addEdge(edgeOnG);
                            uf.unite(fromIndex, toIndex);
                        }
                    }
                }
            }
        }
*/
        System.out.printf("final2 took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
        return this;
    }

    private Graph finalOpr3(Graph g){

        long time = System.nanoTime();
        UnionFind uf = createUnionFind();

        //Gのすべての辺からCに追加できるか試す
        for(Edge gEdge: g.allEdges()) {
            if(!existsEdge(gEdge.from, gEdge.to)){
                int fromI = getVertexIndex(gEdge.from);
                int toI = getVertexIndex(gEdge.to);
                uf.unite(fromI, toI);
                addEdge(gEdge);
            }
            if(uf.isConnected()) break;
        }

        System.out.printf("final3 took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
        return this;
    }

    private String toKey(int from, int to) {
        return "_" + from + "_" + to;
    }

    static private String _toKey(int from, int to) {
        return "_" + from + "_" + to;
    }

    //どの操作でもgraphsを作っているのが冗長
    public Graph opr5(Graph g) {
        long time = System.nanoTime();
        int size = vertices.size();
        if(size < 3 || size > 5) return this;
        LinkedList<Graph> graphs = graphs();
        for(Graph graph: graphs){
            if(graph.isCycle()) continue;
            Vertex[] endpoints = new Vertex[2];
            LinkedList<Vertex> ivs = new LinkedList<>();
            int i = 0;
            for(Vertex vertex: graph.vertices){
                if(vertex.getEdges().size() == 1) {
                    endpoints[i] = vertex;
                    i++;
                }
            }

            int end1 = endpoints[0].getIndex();
            int end2 = endpoints[1].getIndex(); //大小は関係ないはず


            for(int j = 0;j < size;j++){
                //
                if(j == getVertexIndex(end1) || j == getVertexIndex(end2)) continue;
                boolean[] flag = new boolean[size];
                flag[j] = true;
                LinkedList<Edge> path = opr5DFS(g,flag,j,size - 1);
                if(path.isEmpty()) continue;
                for(Vertex v: vertices){
                    v.getEdges().clear();
                }

                for(Edge edge: path){
                    addEdge(edge);
                }
            }
        }
        System.out.printf("opr5 took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
        return this;
    }

    /*
    ivからパスを作る
    prevに初期値を渡す必要があるのでivの数だけ呼び出す必要がある
    呼び出す前にflagの自分の位置をtrueにしないとサイクルになってしまう

    edge + opr5DFS

     */
    public LinkedList<Edge> opr5DFS(Graph g, boolean flag[],int prev, int count){
        if(count == 0) return new LinkedList<>();
        LinkedList<Edge> ret = new LinkedList<>();
        LinkedList<Edge> tmp;
        for(int i = 0;i < getVertices().size();i++){
            if(flag[i]) continue;
            Vertex v = getVertices().get(i);
            if(g.existsEdge(prev, v.getIndex())){
                //Gにおいて辺があるので追加
                Edge edge = new Edge(prev, v.getIndex());
                flag[i] = true;
                tmp = opr5DFS(g, flag, v.getIndex(), count - 1);
                flag[i] = false;
                if(tmp.isEmpty()) continue;
                ret.add(edge);
                ret.addAll(tmp);
                return ret;
            }
        }
        return ret;
    }

    public Graph opr6(Graph g){
        long time = System.nanoTime();
        LinkedList<Graph> paths = new LinkedList<>();
        LinkedList<Graph> cycles = new LinkedList<>();
        for(Graph graph: graphs()){
            if(graph.isCycle()) {
                cycles.add(graph);
            } else {
                paths.add(graph);
            }
        }

        for(Graph path: paths){
            for(Vertex vertex: path.vertices){
                if(vertex.getEdges().size() > 1) continue;
                for(Graph cycle: cycles){
                    for(Vertex cv: cycle.vertices){
                        if(!g.existsEdge(vertex.getIndex(), cv.getIndex())) continue;
                        Edge deleted =  g.getVertex(cv.getIndex()).getEdges().getFirst();
                        deleteEdge(deleted);
                        addEdge(new Edge(vertex.getIndex(), cv.getIndex()));
                    }
                }
            }
        }
        System.out.printf("opr6 took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
        return this;
    }

    public Graph opr7(Graph g){
        long time = System.nanoTime();
        LinkedList<Graph> paths = new LinkedList<>();
        LinkedList<LinkedList<Vertex>> endpointsList = new LinkedList<>();
        LinkedList<LinkedList<Vertex>> ivsList = new LinkedList<>();
        LinkedList<Graph> graphs = graphs();
        //[[1, 2], [5,7]]みたいな
        for(Graph graph: graphs){

            //isPathに変えた方がいいかも
            //PCCなので大丈夫？

            //長さ2のパスはivsが空なので注意
            if(!graph.isCycle()) {
                LinkedList endpoints = new LinkedList();
                LinkedList ivs = new LinkedList();
                //System.out.println("path:");
                for(Vertex vertex: graph.vertices){
                    if(vertex.getEdges().size() == 1) {

                        endpoints.add(vertex);
                        //System.out.println(vertex);
                    } else {
                        ivs.add(vertex);
                    }
                }
                endpointsList.add(endpoints);
                ivsList.add(ivs);
                paths.add(graph);
            }
        }
        System.out.println(endpointsList.size() + ", " + ivsList.size());
        for (int i = 0; i < endpointsList.size(); i++) {
            LinkedList<Vertex> endpoints = endpointsList.get(i);
       // for(LinkedList<Vertex> endpoints: endpointsList){
            int p1Size = pathSize(endpoints.getFirst().getIndex());
            for(Vertex endpoint: endpoints){
                for(int j = 0; j < ivsList.size();j++){
                    //長さ2のパスはivsが空なので注意
                    LinkedList<Vertex> ivs = ivsList.get(j);
                    if(ivs.isEmpty()) continue;
                    if(i == j) continue;//同じパスとは比較しない
                    int p2Size = pathSize(ivs.getFirst().getIndex());
                    for(Vertex iv: ivs){

                        if(!g.existsEdge(endpoint.getIndex(), iv.getIndex())) continue;
                        //辺を付け替えて長さを見る(2回行う)
                        addEdge(new Edge(endpoint.getIndex(), iv.getIndex()));
                        Edge del1 = iv.getEdges().getFirst();
                        Edge del2 = iv.getEdges().getLast();
                        deleteEdge(del1);

                        int q1Size = pathSize(endpointsList.get(j).getFirst().getIndex());
                        int q2Size = pathSize(endpointsList.get(j).getLast().getIndex());
                        //良くなった場合変更を適用し次へ
                        if(Math.max(q1Size, q2Size) > Math.max(p1Size, p2Size)) continue;

                        //2回目
                        addEdge(del1);
                        deleteEdge(del2);

                        q1Size = pathSize(endpointsList.get(j).getFirst().getIndex());
                        q2Size = pathSize(endpointsList.get(j).getLast().getIndex());

                        if(Math.max(q1Size, q2Size) > Math.max(p1Size, p2Size)) continue;
                        //ダメなら元に戻す
                        addEdge(del2);
                    }
                }
            }
        }

        System.out.printf("opr7 took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
        return this;
    }

    public int pathSize(int current){
        int[] count = {0};
        boolean[] isVisited = new boolean[vertices.size()];
        int currentI = getVertexIndex(current);
        return pathSizeDFS(currentI,count,isVisited);
    }

    public int pathSizeDFS(int current, int[] count, boolean[] isVisited){
        for(Edge edge: vertices.get(current).getEdges()) {
            int toI = getVertexIndex(edge.to);
            if(isVisited[toI]) continue;
            isVisited[toI] = true;
            count[0]++;
            pathSizeDFS(toI,count,isVisited);
        }
        return count[0];
    }

    public int numberOfInternalVertices(){
        int i = 0;
        for(Vertex vertex: vertices){
            if(vertex.getEdges().size() > 1) i++;
        }
        return i;
    }

    boolean isConnected(){
        int vSize = vertices.size();
        int[] count = {1};
        boolean[] isVisited = new boolean[vSize];
        isVisited[0] = true;

        return isConnectedDFS(0,count , vSize, isVisited);
    }

    boolean isConnectedDFS(int current,int[] count, int vSize, boolean isVisited[]){
        for(Edge edge: vertices.get(current).getEdges()) {
            int toI = getVertexIndex(edge.to);
            if(isVisited[toI]) continue;
            isVisited[toI] = true;
            count[0]++;
            isConnectedDFS(toI, count, vSize, isVisited);
        }
        return count[0] == vSize;
    }

    private int cycles(){
        int count = 0;
        for(Graph graph: graphs()){
            if(graph.isCycle()) count++;
        }
        return count;
    }

    public Graph opr8() {
        if(vertices.size() <= 3) return this;
        long time = System.nanoTime();
        LinkedList<Vertex> leaves = new LinkedList<Vertex>();
        int lSize = 0;

        for(Vertex vertex: getVertices()) {
            if(vertex.getEdges().size() == 2) {
                leaves.add(vertex);
                lSize++;
            }
        }

        for (int i = 0; i < lSize - 1; i++) {
            for (int j = i + 1; j < lSize; j++) {
                Vertex v3 = leaves.get(i);
                Vertex v4 = leaves.get(j);
                Vertex v3low, v3high;
                LinkedList<Edge> v3Edges = v3.getEdges();
                if(v3Edges.getFirst().to < v3Edges.getLast().to){
                    v3low = getVertex(v3Edges.getFirst().to);
                    v3high = getVertex(v3Edges.getLast().to);
                } else {
                    v3low = getVertex(v3Edges.getLast().to);
                    v3high = getVertex(v3Edges.getFirst().to);
                }

                Vertex v4low, v4high;
                LinkedList<Edge> v4Edges = v4.getEdges();
                if(v4Edges.getFirst().to < v4Edges.getLast().to){
                    v4low = getVertex(v4Edges.getFirst().to);
                    v4high = getVertex(v4Edges.getLast().to);
                } else {
                    v4low = getVertex(v4Edges.getLast().to);
                    v4high = getVertex(v4Edges.getFirst().to);
                }

                if(v3low.getIndex() == v4low.getIndex() && v3high.getIndex() == v4high.getIndex()) {
                    Vertex v1 = v3low;
                    Vertex v2 = v3high;

                    //v1を消して
                }

            }
        }


        System.out.printf("opr8 took %.3f ms to run%n", (System.nanoTime()-time)/1e6);
        return this;
    }
}
