package org.bgi.flexlab.gaea.tools.haplotypecaller.utils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public final class CompressedDataList<T>  implements Iterable<T> {
    protected Map<T,Integer> valueCounts = new HashMap<>();

    public Map<T,Integer> getValueCounts(){
        return valueCounts;
    }

    public boolean isEmpty(){
        return valueCounts.isEmpty();
    }

    @Override
    public Iterator<T> iterator(){
        Iterator<T> it = new Iterator<T>() {
            private Iterator<T> keySetIterator = valueCounts.keySet().iterator();
            private T currentKey = valueCounts.isEmpty() ? null : keySetIterator.next();
            private int currentValueIndex = 0;
            private int currentValueSize = valueCounts.isEmpty() ? 0 : valueCounts.get(currentKey);

            @Override
            public boolean hasNext() {
                return !valueCounts.isEmpty() && (keySetIterator.hasNext() || currentValueIndex < currentValueSize);
            }

            @Override
            public T next() {
                final T retKey = currentKey;
                currentValueIndex++;
                if(currentValueIndex==currentValueSize){
                    if(keySetIterator.hasNext()) {
                        currentKey = keySetIterator.next();
                        currentValueIndex = 0;
                        currentValueSize = valueCounts.get(currentKey);
                    }
                }
                return retKey;
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
        return it;
    }

    @Override
    public String toString(){
        String str = "";
        Object[] keys = valueCounts.keySet().toArray();
        Arrays.sort(keys);
        for (Object i: keys){
            if(!str.isEmpty()) {
                str += ",";
            }
            str+=(i+","+valueCounts.get(i));
        }
        return str;
    }

    public void add(final T val){
        add(val, 1);
    }

    public void add(final T val, final int count){
        if(valueCounts.containsKey(val)){
            valueCounts.put(val, valueCounts.get(val)+count);
        } else {
            valueCounts.put(val, count);
        }

    }

    public void add(final CompressedDataList<T> obj){
        for(final Map.Entry<T, Integer> pair : obj.valueCounts.entrySet()){
            this.add(pair.getKey(), pair.getValue());
        }
    }

}