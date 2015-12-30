package CodingLayer;

public class Coding{
    private int nv;
    private int kv;
    private int wv;
    private int fv;

    private boolean initialized = false;

    public Coding(){
        ;
    }

    public int initialize(int n,int k,int w){
        nv=n;
        kv=k;
        wv=w;
        initialized=true;
        System.out.println("Initialized");
        return 0;
    }
}
