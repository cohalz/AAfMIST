import java.util.ArrayList;
import java.util.LinkedList;
import java.util.TreeSet;

public class UnionFind {
    int[] par;

    UnionFind(int n) {
        par = new int[n];
        for (int i = 0; i < n; i++) par[i] = i;
    }

    public int find(int x) {
        if (par[x] == x) return x;
        return par[x] = find(par[x]);
    }

    public Boolean same(int x, int y) {
        return find(x) == find(y);
    }

    public void unite(int x, int y) {
        if (find(x) == find(y)) return;
        par[find(x)] = find(y);
    }

    public LinkedList<TreeSet<Integer>> getSets() {
        LinkedList<TreeSet<Integer>> lst = new LinkedList<>();
        TreeSet<Integer> tmpSet = new TreeSet<>();
        tmpSet.add(0);
        lst.add(tmpSet);
        for (int i = 1; i < par.length; i++) {
            boolean flag = false;
            for (TreeSet<Integer> set: lst) {
                int e = set.first();
                if(same(i,e)) {
                    flag = true;
                    tmpSet = set;
                }
            }
            if(flag){
                lst.remove(tmpSet);
                tmpSet.add(i);
                lst.add(tmpSet);
            } else {
                tmpSet = new TreeSet<>();
                tmpSet.add(i);
                lst.add(tmpSet);
            }
        }
        return lst;
    }

    public boolean isConnected() {
        for(int i = 0;i < par.length - 1;i++){
            for(int j = i + 1; j < par.length;j++){
                if(!same(i,j)) return false;
            }
        }
        return true;
    }
}