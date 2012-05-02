#ifndef _GOODMAN_KRUSKAL_GAMMA_H_
#define _GOODMAN_KRUSKAL_GAMMA_H_

#include <stdexcept>
#include <iterator>
#include <set>
#include <vector>
#include <algorithm>
#include <stdint.h>

// efficient implementation of Goodman-Kruskal gamma
// 32-bit count_t will easily overflow!

template<typename T, typename count_t = uint64_t>
struct cum_sum_tree{
  // binary tree which maps Y -> n_Y and keeps the cumulative sum of all n_Y with Y_i < Y
  // incomplete, only functions needed for this algorithm.
  // use init with randomly shuffled Y-values to get a more balanced tree.

  cum_sum_tree(): root(NULL) {}
  ~cum_sum_tree(){ delete root; }

  struct node{
    node(T Y_): Y(Y_), n_Y(0), subtree_sum(0), left(NULL), right(NULL) {}
    ~node(){ delete(left); delete(right); }
    
    T Y;
    count_t n_Y;
    count_t subtree_sum;		// sum of n_Y in subtree

    //  private:
    node* left;			// smaller Y
    node* right;		// larger Y
  };


  void init(T Y);
  void add(T Y, count_t increment, count_t& sum_before, count_t& n_equal);

private:
  node* root;
};


template<typename T, typename count_t>
void cum_sum_tree<T, count_t>::init(T Y){
  // insert a node with n_Y = 0

  // search for insertion location
  node** current = &root;
  // node** insert_here = &root;
  while(*current != NULL){
    if(Y < (*current)->Y)
      current = &(*current)->left;
    else if(Y > (*current)->Y)
      current = &(*current)->right;
    else			// ==
      return;			// do nothing, Y is already in the tree
  }

  // insert
  *current = new node(Y);
}


template<typename T, typename count_t>
void cum_sum_tree<T, count_t>::add(T Y, count_t increment, count_t& sum_before, count_t& n_Y){
  // search and update subtree_sum for each encountered node

  sum_before = 0;
  
  node* current = root;
  if(current == 0)
    throw std::runtime_error("Y is not in the cum_sum_tree.");
    
  while(current->Y != Y){
    // add sum of n_Y from subtree with smaller Y
    // if(current->left)
    //   sum_before += current->left->subtree_sum;
    
    current->subtree_sum += increment;
    
    if(Y < current->Y){
      current = current->left;
    }else if(Y > current->Y){
      sum_before += current->n_Y; // because current->Y < Y
      if(current->left)
	sum_before += current->left->subtree_sum;
      current = current->right;
    }

    if(current == 0)
      throw std::runtime_error("Y is not in the cum_sum_tree.");
  }

  current->n_Y += increment;
  current->subtree_sum += increment;

  if(current->left)
    sum_before += current->left->subtree_sum;
  n_Y = current->n_Y;
}


template<typename forward_iterator, typename random_access_iterator>
void secondary_sort(forward_iterator X_begin, forward_iterator X_end, random_access_iterator Y_begin, random_access_iterator Y_end){
  // sort Y for matching X, the pairs X[i], Y[i] will be sorted by X first and then by Y.
  // FIXME: Y_end not necessary, clean up arguments.
  
  size_t i_begin, i_end;
  i_begin = 0; i_end = 1;
  forward_iterator X_previous = X_begin;
  ++X_begin;


  for(forward_iterator i = X_begin; i != X_end; ++i){
    // sort if different or end has been reached
    if(*i != *X_previous){
      // if is wrong:
      // if(i_begin - i_end >= 2)	// sort if necessary
      // if correct:
      if(i_end - i_begin >= 2)
	std::sort(Y_begin + i_begin, Y_begin + i_end);
      i_begin = i_end;
    }

    // final iteration: sort end
    forward_iterator i_2 = i; ++i_2; // tmp iterator
    if(i_2 == X_end)
      std::sort(Y_begin + i_begin, Y_begin + i_end + 1);
    
    ++i_end;
    ++X_previous;
  }
}



template<typename bidirectional_iterator, typename random_access_iterator, typename count_t>
void concordance_count(bidirectional_iterator X_begin, bidirectional_iterator X_end, random_access_iterator Y_begin, random_access_iterator Y_end,
		       count_t& concordant, count_t& discordant, count_t& extraX, count_t& extraY, count_t& spare){
  // For Kendall's tau and goodman-kruskal gamma. Correlation for the (X_i, Y_i) pairs
  // use secondary_sort(...) for the required sorting.
  // X_begin, X_end: sorted range
  // Y_begin, Y_end: sorted for equal X values
  // concordant: X_i > X_j, Y_i > Y_j (or <)
  // discordant: X_i > X_j, Y_i < Y_j (or <, >)
  // extraX:     X_i != X_j, Y_i == Y_j
  // extraY:     X_i == X_j, Y_i != Y_j
  // spare:      X_i == X_j, Y_i == Y_j

  // This is based on algorithm SD from
  // Christensen - Fast algrorithms for the calculation of Kendall's tau, Computational Statistics (2005) 20:51-62

  using namespace std;

  // get unique elements from Y_begin, Y_end by insertion into std::set
  typedef typename iterator_traits<random_access_iterator>::value_type Y_value_type;
  set<Y_value_type> Y_set;
  for(random_access_iterator i = Y_begin; i != Y_end; ++i)
    Y_set.insert(*i);

  // put unique elements into a vector and shuffle randomly for a more balanced tree
  vector<Y_value_type> Y_unique(Y_set.begin(), Y_set.end());
  random_shuffle(Y_unique.begin(), Y_unique.end());

  // init tree
  cum_sum_tree<Y_value_type> cst;
  for(typename vector<Y_value_type>::iterator i = Y_unique.begin(); i != Y_unique.end(); ++i)
    cst.init(*i);
  
  // the counting
  concordant = discordant = extraX = extraY = spare= 0;

  count_t spare_inc = 0;
  count_t extraY_inc = 0;
  // count_t pos = 0;
  // different from X/Y_begin for first iteration of loop:
  typename iterator_traits<bidirectional_iterator>::value_type X_previous = *X_begin + 1;
  typename iterator_traits<random_access_iterator>::value_type Y_previous = *Y_begin + 1;
  random_access_iterator Y_i = Y_begin;
  size_t i = 0;
  
  for(bidirectional_iterator X_i = X_begin; X_i != X_end; ++X_i, ++Y_i, ++i){
    if(*X_i != X_previous){
      extraY_inc = 0;
      spare_inc = 1;
    }else{
      if(*Y_i == Y_previous){
	++spare_inc;
      }else{
	extraY_inc += spare_inc;
	spare_inc = 1;
      }
    }

    count_t sum_before, n_equal;
    cst.add(*Y_i, 1, sum_before, n_equal);

    count_t concordant_inc = sum_before - extraY_inc;
    count_t extraX_inc = n_equal - spare_inc;
    count_t discordant_inc = i - (concordant_inc + extraX_inc + extraY_inc + spare_inc) + 1;

    extraX += extraX_inc;
    extraY += extraY_inc;
    concordant += concordant_inc;
    discordant += discordant_inc;
    spare += spare_inc;

    X_previous = *X_i;
    Y_previous = *Y_i;
  }
}


template<typename count_t>
double kendall_tau(count_t concordant, count_t discordant, count_t extraX, count_t extraY){
  return  ((double) concordant - discordant) / (sqrt(concordant + discordant + extraX) * sqrt(concordant + discordant + extraY));
}

template<typename count_t>
double goodman_kruskal_gamma(count_t concordant, count_t discordant){
  return  ((double) concordant - discordant) / (concordant + discordant);
}

// double kendall_tau_2(size_t concordant, size_t discordant, size_t extraX, size_t extraY){
//   // older definition which does not take ties into account
//   return ((double) concordant - discordant) / (concordant + discordant + extraX + extraY);
// }


#endif /* _GOODMAN_KRUSKAL_GAMMA_H_ */
