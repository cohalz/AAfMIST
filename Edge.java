public class Edge implements Comparable<Edge> {
    public int from;
    public int to;

    public Edge(int from, int to){
        this.from = from;
        this.to = to;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Edge edge = (Edge) o;

        if (from != edge.from) return false;
        return to == edge.to;
    }

    @Override
    public int hashCode() {
        int result = from;
        result = 31 * result + to;
        return result;
    }

    public Edge swapEdge(){
        return new Edge(to, from);
    }

    @Override
    public String toString() {
        return "(" + from + " -> " + to + ")";
    }

    public int compareTo(Edge edge){
        int a = from - edge.from;
        if(a != 0) return a;
        else return to - edge.to;

    }
}
