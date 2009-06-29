open Lacaml.Impl.D

open Interfaces

(* Hyper parameter optimization with the GNU Scientific Library *)
module Gsl = struct

  (* Ordinary hyper parameter optimization *)

  module Make (Spec : Sigs.Deriv) = struct
    open Spec

    let train ?kernel ?sigma2 ?inducing ?n_rand_inducing ~inputs ~targets () =
      let kernel =
        match kernel with
        | None -> Eval.Inputs.create_default_kernel inputs
        | Some kernel -> kernel
      in
      let sigma2 =
        match sigma2 with
        | None -> Vec.sqr_nrm2 targets /. float (Vec.dim targets)
        | Some sigma2 -> max sigma2 min_float
      in
      let eval_inducing_prepared =
        match inducing with
        | None ->
            let n_inducing =
              let n_inputs = Eval.Spec.Inputs.get_n_inputs inputs in
              match n_rand_inducing with
              | None -> min (n_inputs / 10) 1000
              | Some n_rand_inducing -> max (min n_inputs n_rand_inducing) 0
            in
            Eval.Inducing.Prepared.choose_n_random_inputs
              kernel ~n_inducing inputs
        | Some inducing -> Eval.Inducing.Prepared.calc inducing
      in
      let eval_inputs_prepared =
        Eval.Inputs.Prepared.calc eval_inducing_prepared inputs
      in
      let deriv_inducing_prepared =
        Deriv.Inducing.Prepared.calc eval_inducing_prepared
      in
      let deriv_inputs_prepared =
        Deriv.Inputs.Prepared.calc deriv_inducing_prepared eval_inputs_prepared
      in
      let hyper_vars, hyper_vals = Deriv.Spec.Hyper.extract kernel in
      let n_hypers = Array.length hyper_vars in
      let n_gsl_hypers = 1 + n_hypers in
      let gsl_hypers = Gsl_vector.create n_gsl_hypers in
      gsl_hypers.{0} <- log sigma2;
      for i = 1 to n_hypers do gsl_hypers.{i} <- hyper_vals.{i} done;
      let module Gd = Gsl_multimin.Deriv in
      let update_hypers ~gsl_hypers =
        let sigma2 = exp gsl_hypers.{0} in
        let hyper_vals = Vec.create n_hypers in
        for i = 1 to n_hypers do hyper_vals.{i} <- gsl_hypers.{i} done;
        sigma2, Deriv.Spec.Hyper.update kernel hyper_vals
      in
      let multim_f ~x:gsl_hypers =
        let sigma2, kernel = update_hypers ~gsl_hypers in
        let eval_inducing =
          Eval.Inducing.calc kernel eval_inducing_prepared
        in
        let eval_inputs =
          Eval.Inputs.calc eval_inducing eval_inputs_prepared
        in
        let model = Eval.Model.calc eval_inputs ~sigma2 in
        let trained = Eval.Trained.calc model ~targets in
        let log_evidence = Eval.Trained.calc_log_evidence trained in
        -. log_evidence
      in
      let multim_dcommon ~x:gsl_hypers ~g:gradient =
        let sigma2, kernel = update_hypers ~gsl_hypers in
        let deriv_inducing =
          Deriv.Inducing.calc kernel deriv_inducing_prepared
        in
        let deriv_inputs =
          Deriv.Inputs.calc deriv_inducing deriv_inputs_prepared
        in
        let dmodel = Deriv.Model.calc ~sigma2 deriv_inputs in
        let trained = Deriv.Trained.calc dmodel ~targets in
        let dlog_evidence_dsigma2 =
          Deriv.Trained.calc_log_evidence_sigma2 trained
        in
        gradient.{0} <- -. dlog_evidence_dsigma2 *. sigma2;
        let hyper_t = Deriv.Trained.prepare_hyper trained in
        for i = 1 to n_hypers do
          gradient.{i} <-
            -. Deriv.Trained.calc_log_evidence hyper_t hyper_vars.(i - 1)
        done;
        trained
      in
      let multim_df ~x ~g = ignore (multim_dcommon ~x ~g) in
      let multim_fdf ~x ~g =
        let trained = multim_dcommon ~x ~g in
        let log_evidence =
          Eval.Trained.calc_log_evidence (Deriv.Trained.calc_eval trained)
        in
        -. log_evidence
      in
      let multim_fun_fdf =
        {
          Gsl_fun.
          multim_f = multim_f;
          multim_df = multim_df;
          multim_fdf = multim_fdf;
        }
      in
      let mumin =
        Gd.make Gd.VECTOR_BFGS2 n_gsl_hypers
          multim_fun_fdf ~x:gsl_hypers ~step:1e-1 ~tol:1e-1
      in
      let rec loop last_log_evidence =
        let neg_log_likelihood = Gd.minimum ~x:gsl_hypers mumin in
        let log_evidence = -. neg_log_likelihood in
        let diff = abs_float (1. -. (log_evidence /. last_log_evidence)) in
        if diff < 0.001 then
          let sigma2, kernel = update_hypers ~gsl_hypers in
          let eval_inducing =
            Eval.Inducing.calc kernel eval_inducing_prepared
          in
          let eval_inputs =
            Eval.Inputs.calc eval_inducing eval_inputs_prepared
          in
          let model = Eval.Model.calc eval_inputs ~sigma2 in
          Eval.Trained.calc model ~targets
        else begin
          Gd.iterate mumin;
          loop log_evidence
        end
      in
      loop neg_infinity
  end


  (* SPGP *)

  module Make_SPGP
    (Deriv : Sigs.Deriv)
    (Spec : Specs.SPGP
      with module Eval = Deriv.Eval.Spec
      with module Deriv = Deriv.Deriv.Spec) =
  struct
    module Eval = Deriv.Eval
    module Deriv = Deriv

    open Deriv

    module SPGP = struct
      module Spec = Spec

      let train ?kernel ?sigma2 ?inducing ?n_rand_inducing ~inputs ~targets () =
        let kernel =
          match kernel with
          | None -> Eval.Inputs.create_default_kernel inputs
          | Some kernel -> kernel
        in
        let sigma2 =
          match sigma2 with
          | None -> Vec.sqr_nrm2 targets /. float (Vec.dim targets)
          | Some sigma2 -> max sigma2 min_float
        in
        let eval_inducing_prepared =
          match inducing with
          | None ->
              let n_inducing =
                let n_inputs = Eval.Spec.Inputs.get_n_inputs inputs in
                match n_rand_inducing with
                | None -> min (n_inputs / 10) 1000
                | Some n_rand_inducing -> max (min n_inputs n_rand_inducing) 0
              in
              Eval.Inducing.Prepared.choose_n_random_inputs
                kernel ~n_inducing inputs
          | Some inducing -> Eval.Inducing.Prepared.calc inducing
        in
        let hyper_vars, hyper_vals = Deriv.Spec.Hyper.extract kernel in
        let n_hypers = Array.length hyper_vars in
        let inducing_points =
          Eval.Inducing.Prepared.get_points eval_inducing_prepared
        in
        let inducing_hypers, inducing_vals = 
          Spec.Inducing_hypers.extract inducing_points
        in
        let n_inducing_hypers = Array.length inducing_hypers in
        let n_gsl_hypers = 1 + n_hypers + n_inducing_hypers in
        let gsl_hypers = Gsl_vector.create n_gsl_hypers in
        gsl_hypers.{0} <- log sigma2;
        for i = 1 to n_hypers do gsl_hypers.{i} <- hyper_vals.{i} done;
        for i = 1 to n_inducing_hypers do
          gsl_hypers.{n_hypers + i} <- inducing_vals.{i}
        done;
        let module Gd = Gsl_multimin.Deriv in
        let update_hypers ~gsl_hypers =
          let sigma2 = exp gsl_hypers.{0} in
          let hyper_vals = Vec.create n_hypers in
          let inducing_vals = Vec.create n_inducing_hypers in
          for i = 1 to n_hypers do hyper_vals.{i} <- gsl_hypers.{i} done;
          for i = 1 to n_inducing_hypers do
            inducing_vals.{i} <- gsl_hypers.{n_hypers + i}
          done;
          let kernel = Deriv.Spec.Hyper.update kernel hyper_vals in
          let inducing =
            Spec.Inducing_hypers.update inducing_points inducing_vals
          in
          let eval_inducing_prepared = Eval.Inducing.Prepared.calc inducing in
          let eval_inputs_prepared =
            Eval.Inputs.Prepared.calc eval_inducing_prepared inputs
          in
          sigma2, kernel, eval_inducing_prepared, eval_inputs_prepared
        in
        let calc_f_trained gsl_hypers =
          let sigma2, kernel, eval_inducing_prepared, eval_inputs_prepared =
            update_hypers ~gsl_hypers
          in
          let eval_inducing =
            Eval.Inducing.calc kernel eval_inducing_prepared
          in
          let eval_inputs =
            Eval.Inputs.calc eval_inducing eval_inputs_prepared
          in
          let model = Eval.Model.calc eval_inputs ~sigma2 in
          Eval.Trained.calc model ~targets
        in
        let multim_f ~x:gsl_hypers =
          let trained = calc_f_trained gsl_hypers in
          let log_evidence = Eval.Trained.calc_log_evidence trained in
          -. log_evidence
        in
        let multim_dcommon ~x:gsl_hypers ~g:gradient =
          let sigma2, kernel, eval_inducing_prepared, eval_inputs_prepared =
            update_hypers ~gsl_hypers
          in
          let deriv_inducing_prepared =
            Deriv.Inducing.Prepared.calc eval_inducing_prepared
          in
          let deriv_inputs_prepared =
            Deriv.Inputs.Prepared.calc deriv_inducing_prepared
            eval_inputs_prepared
          in
          let deriv_inducing =
            Deriv.Inducing.calc kernel deriv_inducing_prepared
          in
          let deriv_inputs =
            Deriv.Inputs.calc deriv_inducing deriv_inputs_prepared
          in
          let dmodel = Deriv.Model.calc ~sigma2 deriv_inputs in
          let trained = Deriv.Trained.calc dmodel ~targets in
          let dlog_evidence_dsigma2 =
            Deriv.Trained.calc_log_evidence_sigma2 trained
          in
          gradient.{0} <- -. dlog_evidence_dsigma2 *. sigma2;
          let hyper_t = Deriv.Trained.prepare_hyper trained in
          for i = 1 to n_hypers do
            gradient.{i} <-
              -. Deriv.Trained.calc_log_evidence hyper_t hyper_vars.(i - 1)
          done;
          for i = 1 to n_inducing_hypers do
            gradient.{n_hypers + i} <-
              -. Deriv.Trained.calc_log_evidence hyper_t inducing_hypers.(i - 1)
          done;
          trained
        in
        let multim_df ~x ~g = ignore (multim_dcommon ~x ~g) in
        let multim_fdf ~x ~g =
          let trained = multim_dcommon ~x ~g in
          let log_evidence =
            Eval.Trained.calc_log_evidence (Deriv.Trained.calc_eval trained)
          in
          -. log_evidence
        in
        let multim_fun_fdf =
          {
            Gsl_fun.
            multim_f = multim_f;
            multim_df = multim_df;
            multim_fdf = multim_fdf;
          }
        in
        let mumin =
          Gd.make Gd.VECTOR_BFGS2 n_gsl_hypers
            multim_fun_fdf ~x:gsl_hypers ~step:1e-1 ~tol:1e-1
        in
        let rec loop last_log_evidence =
          let neg_log_likelihood = Gd.minimum ~x:gsl_hypers mumin in
          let log_evidence = -. neg_log_likelihood in
          let diff = abs_float (1. -. (log_evidence /. last_log_evidence)) in
          if diff < 0.001 then calc_f_trained gsl_hypers
          else begin
            Gd.iterate mumin;
            loop log_evidence
          end
        in
        loop neg_infinity
    end
  end
end
