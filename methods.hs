forwardEuler :: (Fractional f, Eq f) => (f -> f -> f) -> f -> f -> f -> f -> f
forwardEuler _ _ _ y0 0 = y0
forwardEuler f x0 xn y0 n =
    let h  = (xn - x0) / n
        x1 = x0 + h
        y1 = y0 + (f x0 y0)
    in  forwardEuler f x1 xn y1 (n - 1)

main :: IO ()
main = print (forwardEuler (\_ y -> y) 0 4 1 4 :: Float)
