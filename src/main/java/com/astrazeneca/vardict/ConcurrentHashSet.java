package com.astrazeneca.vardict;

import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

public class ConcurrentHashSet<E> implements Set<E>{


    private Map<E, Boolean> _map = new ConcurrentHashMap<>();

    @Override
    public int size() {
        return _map.size();
    }

    @Override
    public boolean isEmpty() {
        return _map.isEmpty();
    }

    @Override
    public boolean contains(Object e) {
        return _map.containsKey(e);
    }

    @Override
    public Iterator<E> iterator() {
        return _map.keySet().iterator();
    }

    @Override
    public Object[] toArray() {
        return _map.keySet().toArray();
    }

    @Override
    public <T> T[] toArray(T[] a) {
        return _map.keySet().toArray(a);
    }

    @Override
    public boolean add(E e) {
        return _map.put(e, Boolean.TRUE) == null;
    }

    @Override
    public boolean remove(Object o) {
        return _map.remove(o) != null;
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        return _map.keySet().containsAll(c);
    }

    @Override
    public boolean addAll(Collection<? extends E> c) {
        boolean result = false;
        for (E e : c) {
            if (add(e)) {
                result = true;
            }
        }
        return result;
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        return _map.keySet().retainAll(c);
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        return _map.keySet().removeAll(c);
    }

    @Override
    public void clear() {
        _map.clear();
    }

    @Override
    public String toString()
    {
        return _map.keySet().toString();
    }

    @Override
    public int hashCode() {
        return _map.keySet().hashCode();
    }

    @Override
    public boolean equals(Object o)
    {
        return o == this || _map.keySet().equals(o);
    }



}
