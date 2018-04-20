import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

public class Vertex {
    private LinkedList<Edge> edges;

    //小さくする際にインデックスあると便利そう
    private int index;

    public int getIndex() {
        return index;
    }

    @Override
    public String toString() {
        return "index: " + index +", " + edges;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public Vertex(){
        index = -1;
        edges = new LinkedList<>();
    }

    public Vertex(int index){
        this.index = index;
        edges = new LinkedList<>();
    }

    public LinkedList<Edge> getEdges() {
        return edges;
    }

    //反対側の辺追加はされないので注意
    //使う際はかならずfrom < toチェックを使わないように!!!
    public void addEdge(Edge edge) {
        edges.add(edge);
    }

    public List<Edge> filterdEdges(){
        return getEdges().stream()
                .filter(edge -> edge.from < edge.to)
                .collect(Collectors.toList());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Vertex vertex = (Vertex) o;

        if (index != vertex.index) return false;
        return edges != null ? edges.equals(vertex.edges) : vertex.edges == null;
    }

    @Override
    public int hashCode() {
        int result = edges != null ? edges.hashCode() : 0;
        result = 31 * result + index;
        return result;
    }

    public boolean containsEdge(int to){
        for(Edge edge: getEdges()){
            if(edge.to == to) return true;
        }
        return false;
    }

    @Override
    protected Vertex clone() {
        Vertex newVertex = new Vertex(index);
        for(Edge edge: edges) {
            newVertex.addEdge(new Edge(edge.from, edge.to));
        }
        return newVertex;
    }
}
