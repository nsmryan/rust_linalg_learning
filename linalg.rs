
type Vector = Vec<f32>;
type Matrix = Vec<Vector>;


fn solve(a: &Matrix, b: Vector) -> Vector {
    let mut b_matrix = vec!(b);
    let b_col = transpose(&b_matrix);
    println!("b = ");
    print_matrix(&b_matrix);
    println!("bT = ");
    print_matrix(&b_col);
    println!("___");

    let result = matrix_mult(&invert(a), &b_col);

    transpose(&result)[0].clone()
}

fn transpose(m: &Matrix) -> Matrix {
    let mut transposed = Vec::new();

    for j in 0..m[0].len() {
        let mut row = Vec::new();
        for i in 0..m.len() {
            row.push(m[i][j]);
        }
        transposed.push(row);
    }

    transposed
}

fn invert(m: &Matrix) -> Matrix {
    let mut extended = append_matrices(&m, &identity(m.len()));
    let mut reduced = gauss(&extended);
    let split = split_matrix(&reduced, m.len(), 2 * m.len());
    split
}

fn extract_row(m: &Matrix, row: usize) -> Vector { 
    let mut v = Vec::new();

    for i in 0..m[0].len() {
        v.push(m[row][i]);
    }

    v
}

fn extract_col(m: &Matrix, col: usize) -> Vector { 
    let mut v = Vec::new();

    for i in 0..m.len() {
        v.push(m[i][col]);
    }

    v
}


fn vector_dot(v: &Vector, w: &Vector) -> f32 {
    let mut result = 0.0;

    for index in 0..v.len() {
        result += v[index] * w[index];
    }

    result
}

fn row_dot_col(m: &Matrix, n: &Matrix, row_index: usize, col_index: usize) -> f32 {
    let mut result : f32 = 0.0;

    for i in 0..m[0].len() {
        result += m[row_index][i] * n[i][col_index];
    }

    result
}

fn matrix_mult(m: &Matrix, n: &Matrix) -> Matrix {
    let mut matrix = Vec::new();

    for row_index in 0..m.len() {
        let mut row = Vec::new();

        for col_index in 0..n[0].len() {
            let row_vector = extract_row(m, row_index);
            let col_vector = extract_col(n, col_index);
            //println!("row:");
            //print_vector(&row_vector);
            //println!("col:");
            //print_vector(&col_vector);
            let result = vector_dot(&row_vector, &col_vector);
            //println!("result = {}", result);
            row.push(result);
        }
        matrix.push(row);
    }

    matrix
}

fn addmult(matrix: &mut Matrix, src_row: usize, dst_row: usize, multiplier: f32) {
    for index in 0..matrix[0].len() {
        matrix[dst_row][index] -= matrix[src_row][index] * multiplier;
    }
}

fn scale_matrix(matrix: &mut Matrix, scaler: f32) {
    for i in 0..matrix.len() {
        for j in 0..matrix[0].len() {
            matrix[i][j] *= scaler;
        }
    }
}

fn scale_row(matrix: &mut Matrix, row: usize, scaler: f32) {
    for i in 0..matrix[0].len() {
      matrix[row][i] *= scaler;
    }
}

fn swap_rows(matrix: &mut Matrix, src_row: usize, dst_row: usize) {
    for index in 0..matrix[src_row].len() {
        let tmp = matrix[dst_row][index];

        matrix[dst_row][index] = matrix[src_row][index];
        matrix[src_row][index] = tmp;
    }
}

fn gauss(m: &Matrix) -> Matrix {
    let mut matrix = m.clone();

    //println!("starting gaussian elimination:");
    for i in 0..(matrix.len() - 1) {
        // check if the entry is 0. If so, attempt to pivot.
        if matrix[i][i] == 0.0 {
            // look over the remaining rows in case there is a non-zero entry
            // NOTE this could look for the smallest entry for additional numeric
            // stability
            for k in i+1..matrix.len() {
                if matrix[k][i] != 0.0 {
                    swap_rows(&mut matrix, i, k);
                    break;
                }
            }
        }

        for j in i+1..matrix.len() {
            //print_matrix(&matrix);
            //println!("");
            let multiplier: f32 = matrix[j][i] / matrix[i][i];

           addmult(&mut matrix, i, j, multiplier)
        }
    }
    //print_matrix(&matrix);
    //println!("ending gaussian elimination:");

    matrix
}

fn upper_elim(matrix: &mut Matrix) {
    for i in 1..matrix.len() {
        for j in 0..i {
           let multiplier: f32 = matrix[j][i] / matrix[i][i];

           addmult(matrix, i, j, multiplier)
        }
    }
}

fn ones_diag(matrix: &mut Matrix) {
    for i in 0..matrix.len() {
        if matrix[i][i] != 1.0 {
            let scaler = 1.0 / matrix[i][i];
            scale_row(matrix, i, scaler);
        }
    }
}

fn split_matrix(matrix: &Matrix, start_col: usize, end_col: usize) -> Matrix {
    let mut m: Matrix = Vec::new();

    for i in 0..matrix.len() {
        let mut row = Vec::new();
        for j in start_col..end_col {
            row.push(matrix[i][j]);
        }
        m.push(row);
    }

    m
}

fn identity(size: usize) -> Matrix {
    let mut identity = Vec::new();

    for i in 0..size {
        let mut row = Vec::new();
        for j in 0..size {
            if i == j {
                row.push(1.0);
            }
            else {
                row.push(0.0);
            }
        }
        identity.push(row);
    }

    identity
}

fn append_matrices(m: &Matrix, n: &Matrix) -> Matrix { 
    let mut matrix = Vec::new();

    for row_index in 0..m.len() {
        let mut row = Vec::new();
        row.extend(&m[row_index]);
        row.extend(&n[row_index]);
        matrix.push(row);
    }

    matrix
}

fn print_matrix(m: &Matrix) {
    for i in 0..m.len() {
        print_vector(&m[i]);
    }
}

fn print_vector(v: &Vector) {
    for i in 0..v.len() {
        print!("{:4.1} ", v[i]);
    }
    println!("")
}


fn main () {
    let m: Matrix = vec!(vec!(1.0, 2.0, 3.0),
                         vec!(4.0, 5.0, 6.0),
                         vec!(7.0, 8.0, 10.0));
    println!("matrix:");
    print_matrix(&m);

    println!("transposed matrix:");
    let trans = transpose(&m);
    print_matrix(&trans);

    println!("gaussian elimination:");
    let mut m_prime = gauss(&m);
    print_matrix(&m_prime);

    println!("extended:");
    let extended = append_matrices(&m, &identity(3));
    print_matrix(&extended);

    println!("extended and gaussian:");
    let mut extended_prime = gauss(&extended);
    print_matrix(&extended_prime);

    println!("ones on diag");
    ones_diag(&mut extended_prime);
    print_matrix(&extended_prime);

    println!("upper elimination:");
    upper_elim(&mut extended_prime);
    print_matrix(&extended_prime);

    println!("split out inverted matrix:");
    let inverted = split_matrix(&extended_prime, 3, 6);
    print_matrix(&inverted);

    println!("original:");
    print_matrix(&m);

    println!("multiply by original:");
    let iden = matrix_mult(&m, &inverted);
    print_matrix(&iden);

    println!("multiply by original, other:");
    let iden_other = matrix_mult(&inverted, &m);
    print_matrix(&iden_other);

    println!("solve Ax = b");
    let solution = solve(&m, vec!(1.0, 2.0, 3.0));
    print_vector(&solution);
}
