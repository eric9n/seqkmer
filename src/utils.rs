#[derive(Debug, Clone)]
pub enum OptionPair<T> {
    Single(T),
    Pair(T, T),
}

impl<T> OptionPair<T> {
    /// Returns a reference to the single value if it exists.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::OptionPair;
    ///
    /// let single = OptionPair::Single(42);
    /// assert_eq!(single.single(), Some(&42));
    ///
    /// let pair = OptionPair::Pair(1, 2);
    /// assert_eq!(pair.single(), None);
    /// ```
    pub fn single(&self) -> Option<&T> {
        match self {
            OptionPair::Single(value) => Some(value),
            _ => None,
        }
    }

    /// Maps the OptionPair using a provided function.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::OptionPair;
    ///
    /// let single = OptionPair::Single(42);
    /// let mapped_single = single.map(|x| Ok::<_, ()>(x * 2)).unwrap();
    /// assert_eq!(mapped_single, OptionPair::Single(84));
    ///
    /// let pair = OptionPair::Pair(1, 2);
    /// let mapped_pair = pair.map(|x| Ok::<_, ()>(x * 2)).unwrap();
    /// assert_eq!(mapped_pair, OptionPair::Pair(2, 4));
    /// ```
    pub fn map<U, E, F>(&self, mut f: F) -> Result<OptionPair<U>, E>
    where
        F: FnMut(&T) -> Result<U, E>,
    {
        match self {
            OptionPair::Single(t) => f(t).map(OptionPair::Single),
            OptionPair::Pair(t1, t2) => {
                let u1 = f(t1)?;
                let u2 = f(t2)?;
                Ok(OptionPair::Pair(u1, u2))
            }
        }
    }

    /// Reduces the OptionPair to a single value using a provided function.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::OptionPair;
    ///
    /// let single = OptionPair::Single(42);
    /// let reduced_single = single.reduce(0, |acc, &x| acc + x);
    /// assert_eq!(reduced_single, 42);
    ///
    /// let pair = OptionPair::Pair(1, 2);
    /// let reduced_pair = pair.reduce(0, |acc, &x| acc + x);
    /// assert_eq!(reduced_pair, 3);
    /// ```
    pub fn reduce<U, F>(&self, init: U, mut f: F) -> U
    where
        F: FnMut(U, &T) -> U,
    {
        match self {
            OptionPair::Single(t) => f(init, t),
            OptionPair::Pair(t1, t2) => {
                let result = f(init, t1);
                f(result, t2)
            }
        }
    }

    /// Reduces the OptionPair to a string using a provided function and separator.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::OptionPair;
    ///
    /// let single = OptionPair::Single(42);
    /// let reduced_single = single.reduce_str(", ", |&x| x.to_string());
    /// assert_eq!(reduced_single, "42");
    ///
    /// let pair = OptionPair::Pair(1, 2);
    /// let reduced_pair = pair.reduce_str(", ", |&x| x.to_string());
    /// assert_eq!(reduced_pair, "1, 2");
    /// ```
    pub fn reduce_str<F>(&self, sep: &str, mut f: F) -> String
    where
        F: FnMut(&T) -> String,
    {
        self.reduce(String::new(), |acc, t| {
            if acc.is_empty() {
                f(t)
            } else {
                format!("{}{}{}", acc, sep, f(t))
            }
        })
    }

    /// Applies a function to each value in the OptionPair.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::OptionPair;
    ///
    /// let single = OptionPair::Single(42);
    /// let applied_single = single.apply(|&x| x * 2);
    /// assert_eq!(applied_single, OptionPair::Single(84));
    ///
    /// let pair = OptionPair::Pair(1, 2);
    /// let applied_pair = pair.apply(|&x| x * 2);
    /// assert_eq!(applied_pair, OptionPair::Pair(2, 4));
    /// ```
    pub fn apply<U, F>(&self, mut f: F) -> OptionPair<U>
    where
        F: FnMut(&T) -> U,
    {
        match self {
            OptionPair::Single(t) => OptionPair::Single(f(t)),
            OptionPair::Pair(t1, t2) => OptionPair::Pair(f(t1), f(t2)),
        }
    }

    /// Applies a mutable function to each value in the OptionPair.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::OptionPair;
    ///
    /// let mut single = OptionPair::Single(42);
    /// let applied_single = single.apply_mut(|x| *x * 2);
    /// assert_eq!(applied_single, OptionPair::Single(84));
    ///
    /// let mut pair = OptionPair::Pair(1, 2);
    /// let applied_pair = pair.apply_mut(|x| *x * 2);
    /// assert_eq!(applied_pair, OptionPair::Pair(2, 4));
    /// ```
    pub fn apply_mut<U, F>(&mut self, mut f: F) -> OptionPair<U>
    where
        F: FnMut(&mut T) -> U,
    {
        match self {
            OptionPair::Single(t) => OptionPair::Single(f(t)),
            OptionPair::Pair(t1, t2) => OptionPair::Pair(f(t1), f(t2)),
        }
    }
}

impl<T: Clone> OptionPair<T> {
    /// Creates an OptionPair from a slice.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::OptionPair;
    ///
    /// let single = OptionPair::from_slice(&[42]);
    /// assert_eq!(single, OptionPair::Single(42));
    ///
    /// let pair = OptionPair::from_slice(&[1, 2]);
    /// assert_eq!(pair, OptionPair::Pair(1, 2));
    /// ```
    pub fn from_slice(slice: &[T]) -> OptionPair<T> {
        match slice {
            [a, b] => OptionPair::Pair(a.clone(), b.clone()),
            [a] => OptionPair::Single(a.clone()),
            _ => unreachable!(),
        }
    }
}

impl<T> From<(T, Option<T>)> for OptionPair<T> {
    /// Creates an OptionPair from a tuple of (T, Option<T>).
    ///
    /// # Examples
    ///
    /// ```
    /// use seqkmer::OptionPair;
    ///
    /// let single: OptionPair<i32> = (42, None).into();
    /// assert_eq!(single, OptionPair::Single(42));
    ///
    /// let pair: OptionPair<i32> = (1, Some(2)).into();
    /// assert_eq!(pair, OptionPair::Pair(1, 2));
    /// ```
    fn from(tuple: (T, Option<T>)) -> Self {
        match tuple {
            (a, Some(b)) => OptionPair::Pair(a, b),
            (a, None) => OptionPair::Single(a),
        }
    }
}

impl<T: PartialEq> PartialEq for OptionPair<T> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (OptionPair::Single(a), OptionPair::Single(b)) => a == b,
            (OptionPair::Pair(a1, a2), OptionPair::Pair(b1, b2)) => a1 == b1 && a2 == b2,
            _ => false,
        }
    }
}
