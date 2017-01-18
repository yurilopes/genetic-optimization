template <typename Iter>
Iter nextIterator(Iter iter){
    return ++iter;
}

template <typename Iter, typename Cont>
bool isLastIterator(Iter iter, const Cont &cont){
    return (iter != cont.end()) && (nextIterator(iter) == cont.end());
}
