open Lacaml.Impl.D
open Interfaces

(* Hyper parameter optimization with the GNU Scientific Library *)
module Gsl = struct
  exception Optim_exception of exn

  let check_exception seen_exception_ref res =
    if classify_float res = FP_nan then
      match !seen_exception_ref with
      | None -> failwith "Gpr.Optim.Gsl: optimization function returned nan"
      | Some exc -> raise (Optim_exception exc)

  let ignore_report ~iter:_ _ = ()

  let get_n_rand_inducing loc n_inputs = function
    | None -> min (n_inputs / 10) 1000
    | Some n_rand_inducing ->
        if n_rand_inducing < 1 then
          failwith (
            Printf.sprintf
              "Gpr.Optim%s.train: n_rand_inducing (%d) < 1" loc
                n_rand_inducing)
        else if n_rand_inducing > n_inputs then
          failwith (
            Printf.sprintf
              "Gpr.Optim%s.train: n_rand_inducing (%d) > n_inputs (%d)" loc
                n_rand_inducing n_inputs)
        else n_rand_inducing


  (* Ordinary hyper parameter optimization *)

  module Make (Spec : Sigs.Deriv) = struct
    open Spec

    let train
          ?(step = 1e-1) ?(tol = 1e-1) ?(epsabs = 0.1)
          ?(report_trained_model = ignore_report)
          ?(report_gradient_norm = (fun ~iter:_ _ -> ()))
          ?kernel ?sigma2 ?inducing ?n_rand_inducing ~inputs ~targets () =
      let sigma2 =
        match sigma2 with
        | None -> Vec.sqr_nrm2 targets /. float (Vec.dim targets)
        | Some sigma2 -> max sigma2 min_float
      in
      let kernel, inducing =
        match inducing with
        | None ->
            let n_inducing =
              let n_inputs = Eval.Spec.Inputs.get_n_points inputs in
              get_n_rand_inducing "" n_inputs n_rand_inducing
            in
            let kernel =
              match kernel with
              | None -> Eval.Inputs.create_default_kernel ~n_inducing inputs
              | Some kernel -> kernel
            in
            (
              kernel,
              Eval.Inducing.choose_n_random_inputs kernel ~n_inducing inputs
            )
        | Some inducing ->
            match kernel with
            | None ->
                let n_inducing = Eval.Spec.Inducing.get_n_points inducing in
                Eval.Inputs.create_default_kernel ~n_inducing inputs, inducing
            | Some kernel -> kernel, inducing
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
      let seen_exception_ref = ref None in
      let wrap_seen_exception f =
        try f () with exc -> seen_exception_ref := Some exc; raise exc
      in
      let best_model_ref = ref None in
      let get_best_model () =
        match !best_model_ref with
        | None -> assert false  (* impossible *)
        | Some (trained, _) -> trained
      in
      let iter_count = ref 1 in
      let update_best_model trained log_evidence =
        match !best_model_ref with
        | Some (_, old_log_evidence)
          when old_log_evidence >= log_evidence -> ()
        | _ ->
            report_trained_model ~iter:!iter_count trained;
            best_model_ref := Some (trained, log_evidence)
      in
      let multim_f ~x:gsl_hypers =
        let sigma2, kernel = update_hypers ~gsl_hypers in
        let eval_inducing = Eval.Inducing.calc kernel inducing in
        let eval_inputs = Eval.Inputs.calc eval_inducing inputs in
        let model = Eval.Model.calc eval_inputs ~sigma2 in
        let trained = Eval.Trained.calc model ~targets in
        let log_evidence = Eval.Trained.calc_log_evidence trained in
        update_best_model trained log_evidence;
        -. log_evidence
      in
      let multim_f ~x = wrap_seen_exception (fun () -> multim_f ~x) in
      let multim_dcommon ~x:gsl_hypers ~g:gradient =
        let sigma2, kernel = update_hypers ~gsl_hypers in
        let deriv_inducing = Deriv.Inducing.calc kernel inducing in
        let deriv_inputs = Deriv.Inputs.calc deriv_inducing inputs in
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
      let multim_df ~x ~g = wrap_seen_exception (fun () -> multim_df ~x ~g) in
      let multim_fdf ~x ~g =
        let deriv_trained = multim_dcommon ~x ~g in
        let trained = Deriv.Trained.calc_eval deriv_trained in
        let log_evidence = Eval.Trained.calc_log_evidence trained in
        update_best_model trained log_evidence;
        -. log_evidence
      in
      let multim_fdf ~x ~g =
        wrap_seen_exception (fun () -> multim_fdf ~x ~g)
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
          multim_fun_fdf ~x:gsl_hypers ~step ~tol
      in
      let gsl_dhypers = Gsl_vector.create n_gsl_hypers in
      let rec loop () =
        let neg_log_likelihood =
          Gd.minimum ~x:gsl_hypers ~g:gsl_dhypers mumin
        in
        check_exception seen_exception_ref neg_log_likelihood;
        let gnorm = Gsl_blas.nrm2 gsl_dhypers in
        report_gradient_norm ~iter:!iter_count gnorm;
        if gnorm < epsabs then get_best_model ()
        else begin
          incr iter_count;
          Gd.iterate mumin;
          loop ()
        end
      in
      loop ()
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

      module Gsl = struct
        let train
              ?(step = 1e-1) ?(tol = 1e-1) ?(epsabs = 0.1)
              ?(report_trained_model = ignore_report)
              ?(report_gradient_norm = (fun ~iter:_ _ -> ()))
              ?kernel ?sigma2 ?inducing ?n_rand_inducing ~inputs ~targets () =
          let sigma2 =
            match sigma2 with
            | None -> Vec.sqr_nrm2 targets /. float (Vec.dim targets)
            | Some sigma2 -> max sigma2 min_float
          in
          let kernel, inducing =
            match inducing with
            | None ->
                let n_inducing =
                  let n_inputs = Eval.Spec.Inputs.get_n_points inputs in
                  get_n_rand_inducing ".SPGP" n_inputs n_rand_inducing
                in
                let kernel =
                  match kernel with
                  | None -> Eval.Inputs.create_default_kernel ~n_inducing inputs
                  | Some kernel -> kernel
                in
                (
                  kernel,
                  Eval.Inducing.choose_n_random_inputs kernel ~n_inducing inputs
                )
            | Some inducing ->
                match kernel with
                | None ->
                    let n_inducing = Eval.Spec.Inducing.get_n_points inducing in
                    (
                      Eval.Inputs.create_default_kernel ~n_inducing inputs,
                      inducing
                    )
                | Some kernel -> kernel, inducing
          in
          let hyper_vars, hyper_vals = Deriv.Spec.Hyper.extract kernel in
          let n_hypers = Array.length hyper_vars in
          let inducing_hypers, inducing_vals =
            Spec.Inducing_hypers.extract inducing
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
            let inducing = Spec.Inducing_hypers.update inducing inducing_vals in
            sigma2, kernel, inducing
          in
          let seen_exception_ref = ref None in
          let wrap_seen_exception f =
            try f () with exc -> seen_exception_ref := Some exc; raise exc
          in
          let best_model_ref = ref None in
          let get_best_model () =
            match !best_model_ref with
            | None -> assert false  (* impossible *)
            | Some (trained, _) -> trained
          in
          let iter_count = ref 1 in
          let update_best_model trained log_evidence =
            match !best_model_ref with
            | Some (_, old_log_evidence)
              when old_log_evidence >= log_evidence -> ()
            | _ ->
                report_trained_model ~iter:!iter_count trained;
                best_model_ref := Some (trained, log_evidence)
          in
          let calc_f_trained gsl_hypers =
            let sigma2, kernel, inducing = update_hypers ~gsl_hypers in
            let eval_inducing = Eval.Inducing.calc kernel inducing in
            let eval_inputs = Eval.Inputs.calc eval_inducing inputs in
            let model = Eval.Model.calc eval_inputs ~sigma2 in
            Eval.Trained.calc model ~targets
          in
          let multim_f ~x:gsl_hypers =
            let trained = calc_f_trained gsl_hypers in
            let log_evidence = Eval.Trained.calc_log_evidence trained in
            update_best_model trained log_evidence;
            -. log_evidence
          in
          let multim_f ~x = wrap_seen_exception (fun () -> multim_f ~x) in
          let multim_dcommon ~x:gsl_hypers ~g:gradient =
            let sigma2, kernel, inducing =
              update_hypers ~gsl_hypers
            in
            let deriv_inducing = Deriv.Inducing.calc kernel inducing in
            let deriv_inputs = Deriv.Inputs.calc deriv_inducing inputs in
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
                -.
                  Deriv.Trained.calc_log_evidence
                    hyper_t inducing_hypers.(i - 1)
            done;
            trained
          in
          let multim_df ~x ~g = ignore (multim_dcommon ~x ~g) in
          let multim_df ~x ~g =
            wrap_seen_exception (fun () -> multim_df ~x ~g)
          in
          let multim_fdf ~x ~g =
            let deriv_trained = multim_dcommon ~x ~g in
            let trained = Deriv.Trained.calc_eval deriv_trained in
            let log_evidence = Eval.Trained.calc_log_evidence trained in
            update_best_model trained log_evidence;
            -. log_evidence
          in
          let multim_fdf ~x ~g =
            wrap_seen_exception (fun () -> multim_fdf ~x ~g)
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
              multim_fun_fdf ~x:gsl_hypers ~step ~tol
          in
          let gsl_dhypers = Gsl_vector.create n_gsl_hypers in
          let rec loop () =
            let neg_log_likelihood =
              Gd.minimum ~x:gsl_hypers ~g:gsl_dhypers mumin
            in
            check_exception seen_exception_ref neg_log_likelihood;
            let gnorm = Gsl_blas.nrm2 gsl_dhypers in
            report_gradient_norm ~iter:!iter_count gnorm;
            if gnorm < epsabs then get_best_model ()
            else begin
              incr iter_count;
              Gd.iterate mumin;
              loop ()
            end
          in
          loop ()
      end
    end
  end
end
